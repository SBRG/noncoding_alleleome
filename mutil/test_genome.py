from genome import get_BOP27_pos_from_K12_pos, is_overlap, get_promoter_range_from_RegulonDB_df_row 

assert(get_BOP27_pos_from_K12_pos(1)==1)
assert(get_BOP27_pos_from_K12_pos(3560456)==3560455)
assert(get_BOP27_pos_from_K12_pos(4296381)==4296382)
assert(get_BOP27_pos_from_K12_pos(2173363)==-1)
assert(get_BOP27_pos_from_K12_pos(2173364)==-1)

assert(is_overlap((1,10),(8,20))==True)
assert(is_overlap((1,10),(4,5))==True)
assert(is_overlap((1,10),(11,12))==False)

# RegulonDB defines promoters as generally being 80 BP in length.
# The unit test demonstrates that what is returned is 80 BP in lenght
assert(get_promoter_range_from_RegulonDB_df_row([0,0,"forward",70,0,0,0,0,0,0]) == (10,90)) 
assert(get_promoter_range_from_RegulonDB_df_row([0,0,"reverse",70,0,0,0,0,0,0]) == (50,130))

print("DONE")
