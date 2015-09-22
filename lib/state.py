from bm_ops import *
from smp_data_structures import Component
import logging
##
# This module is used to store and process state information
# in the semi-Markov simulation.
#
# Kevin Greenan (kmgreen@cs.ucsc.edu)
#

##
# Print debug messages
#
#logging.basicConfig(level=logging.DEBUG)

##
# This encapsulated system state of the components
#
class State:

	CURR_STATE_OK = "system is operational"
	CURR_STATE_DEGRADED = "system has at least one failure"

	##
	#  Given a list of component IDs construct the data 
	#  structures needed to capture system state.
	#  
	def __init__(self, components=[]):
		self.num_components = len(components)
		self.components = components

		# Keep track of number of component failures
		self.num_failed_comp = 0
		
		# Bit-map of failed components
		self.failed_comp = 0

		# Bit-map of available components
		self.avail_comp = (1 << self.num_components) - 1

		# System state
		self.sys_state = self.CURR_STATE_OK
	
	##
	# Function to copy another State Obj
	#
	def copy(self, state):
		self.num_components = state.num_components
		self.components = state.components[:]
		self.num_failed_comp = state.num_failed_comp
		self.failed_comp = state.failed_comp
		self.avail_comp = state.avail_comp
		self.sys_state = state.sys_state

	##
	# Generic update call for state transitions
	#
	def update_state(self, event_type, args):
		if event_type == Component.EVENT_COMP_FAIL:
			self.fail_component(args[0])
			self.sys_state = self.CURR_STATE_DEGRADED
		if event_type == Component.EVENT_COMP_REPAIR:
			self.repair_component(args[0])
			if self.num_failed_comp == 0:
				self.sys_state = self.CURR_STATE_OK
			
		else:
			return None

		#
		# Return a state summary: num failed components + comp
		#
		return self.get_num_component_fail()
		

	##
	# Set a component as failed 
	#
	def fail_component(self, component_id):
		# Insert into bitmap of component failures
		self.failed_comp = bm_insert(self.failed_comp, component_id)
		self.avail_comp = bm_rm(self.avail_comp, component_id)	
		
		# Incrememnt component failure count
		self.num_failed_comp += 1

		# Logging message
		logging.debug("Component %s has failed" % component_id)

	##
	# "Repair" a component by removing it from the failed component list
	#
	def repair_component(self, component_id):
		# Remove entry from the failed component bitmap
		self.failed_comp = bm_rm(self.failed_comp, component_id)	
		self.avail_comp = bm_insert(self.avail_comp, component_id)

		# Decrement component failure count
		self.num_failed_comp -= 1
		
		# Logging message
		logging.debug("Component %s has been repaired" % component_id)

	##
	# Get number of component failures
	#
	def get_num_component_fail(self):
		return self.num_failed_comp

	##
	# Return a list of component IDs of failed components
	#
	def get_failed_components(self):
		return bm_to_list(self.failed_comp)

	##
	# Return a list of component IDs of available components
	#
	def get_avail_components(self):
		return bm_to_list(self.avail_comp)

	##
	# Return the current system state
	#
	def get_sys_state(self):
		return self.sys_state

