#
# Generate {len(list) \choose k} combinations of the elements of
# list.  Result is put in all_combs
#
def generate_combinations(list, k, all_combs, tmp = []):
    if k == 0:
        all_combs.append(tmp[:])
    else:
        for i in range(len(list)):
            tmp.append(list[i])
            generate_combinations(list[i+1:], k-1, all_combs, tmp)
            tmp.pop()

# all combos of size k
# pass in ndx=0, kset=[], and l_kset=[] on first call. 
def kset_set(set, k, ndx, kset, l_kset):
    if len(kset) == k:
        l_kset.append(kset)
        return l_kset
    if ndx == len(set):
        return l_kset
    for i in range(ndx,len(set)):
        nxt_kset = kset[:]
        nxt_kset.append(set[i])
        l_kset = kset_set(set, k, i+1, nxt_kset, l_kset)
    return l_kset

