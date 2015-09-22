from numpy import random
from mpmath import exp
from mpmath import mpf
from mpmath import power
from mpmath import mp
from mpmath import binomial
from plotting import DataPoints, Plotting

mp.prec += 75

class SectorFailModel:
    def __init__(self, total_num_sectors, sectors_per_region, scrub_interval, sector_fail_prob, request_rate=1, write_ratio=1):
        self.total_num_sectors = total_num_sectors

        # If scrub interval is not given, ignore all of the scrubbing-specific parameters
        if scrub_interval is not None:
            self.scrub_interval = scrub_interval
            self.scrub_rate = mpf(1)/self.scrub_interval
            self.sectors_per_region = sectors_per_region
            self.num_regions = self.total_num_sectors / self.sectors_per_region
            self.alpha = (mpf(1)/request_rate) * sector_fail_prob
            self.disk_scrub_period = (mpf(scrub_interval) / sectors_per_region) * total_num_sectors
        
        self.request_rate = request_rate

        self.sector_fail_prob = sector_fail_prob
        
        self.write_ratio = write_ratio
    
    def prob_of_bad_sector(self, time=0):
        return 0 

class RandomScrubSectorFailModel(SectorFailModel):
    
    def prob_of_bad_sector(self, time=0):
        prob = ((self.request_rate * self.disk_scrub_period) / (1+(self.request_rate * self.disk_scrub_period))) * (self.sector_fail_prob*self.write_ratio)

        return 1 - ((1 - prob)**self.total_num_sectors)
    
class DeterministicScrubSectorFailModel(SectorFailModel):

    def prob_of_bad_sector(self, time=0):
        prob = (1 - ((1-exp(-self.request_rate * self.disk_scrub_period))/(self.request_rate * self.disk_scrub_period)))*(self.sector_fail_prob*self.write_ratio)
        
        return 1 - ((1 - prob)**self.total_num_sectors)
    
class BERSectorFailModel(SectorFailModel):
        
    def prob_of_bad_sector(self, time=0):
        prob = self.sector_fail_prob
        return 1 - ((1 - self.sector_fail_prob)**self.total_num_sectors)
        #return (self.sector_fail_prob*self.total_num_sectors)

class NoScrubSectorFailModel(SectorFailModel):
        
    def prob_of_bad_sector(self, time=0):
        prob = self.write_ratio * self.sector_fail_prob
        return 1 - ((1 - prob)**self.total_num_sectors)
    
def test():
    bytes_in_terabyte = 1099511627776
    scrub_fraction = 1
    sector_size = 512
    num_sectors_per_disk = (585937500) * 1
    num_tb_per_disk = mpf(sector_size * num_sectors_per_disk) / bytes_in_terabyte
    sector_error_rate = 4.096e-11
    sector_writes_per_disk_per_hr = mpf(143750) / 24
    writes_to_sector_per_hour = 0.0045
    low_scrub = 168
    high_scrub = 169
    scrub_incr = 12
    
    print writes_to_sector_per_hour
    
    x = []
    y = []
    dp = DataPoints()
    for i in range(low_scrub, high_scrub, scrub_incr):
        ssfm = RandomScrubSectorFailModel(num_sectors_per_disk, scrub_fraction*num_sectors_per_disk, i, sector_error_rate, writes_to_sector_per_hour, 1)
        #x.append(ssfm.disk_scrub_period)
        x.append(i)
        y.append(ssfm.prob_of_bad_sector(i))
        
    dp.addDataSet(x,y,"Random Scrub Model")
    
    x = []
    y = []
    for i in range(low_scrub, high_scrub, scrub_incr):
        
        dsfm = DeterministicScrubSectorFailModel(num_sectors_per_disk, scrub_fraction*num_sectors_per_disk, i, sector_error_rate, writes_to_sector_per_hour, 1)
        #x.append(ssfm.disk_scrub_period)
        x.append(i)
        y.append(dsfm.prob_of_bad_sector(i))
        print "Scrub", i, dsfm.prob_of_bad_sector(i)
        
    dp.addDataSet(x,y,"Deterministic Scrub Model")
    
    x = []
    y = []
    
    for i in range(low_scrub, high_scrub, scrub_incr):
        berfm = BERSectorFailModel(num_sectors_per_disk, scrub_fraction*num_sectors_per_disk, i, sector_error_rate, writes_to_sector_per_hour, 1)
        #x.append(ssfm.disk_scrub_period)
        x.append(i)
        y.append(berfm.prob_of_bad_sector())
        print "BER", i, berfm.prob_of_bad_sector()
    dp.addDataSet(x,y,"BER Model")
    
    x = []
    y = []
    
    for i in range(low_scrub, high_scrub, scrub_incr):
        nsfm = NoScrubSectorFailModel(num_sectors_per_disk, scrub_fraction*num_sectors_per_disk, i, sector_error_rate, writes_to_sector_per_hour, 1)
        #x.append(ssfm.disk_scrub_period)
        x.append(i)
        y.append(nsfm.prob_of_bad_sector())
    
    dp.addDataSet(x,y,"No Scrub Model")
    
    
    pl = Plotting(dp, "graphs/sector_fail_write_%.1f.%.3f.pdf" % (num_tb_per_disk, scrub_fraction), "Comparison of Sector Failure Models \n for %.1f TB Disk (Region %d Perc.)" % (num_tb_per_disk, 100*scrub_fraction), "Scrub Interval (hours)", "Probability of Sector Failure", type=Plotting.YLOG, legend_loc=2)
    pl.plot()

if __name__ == "__main__":
    test()
