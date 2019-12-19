
# In[1]:

import pickle

last_col_path='./../data/chrX_last_col.txt'
chrX_map_path='./../data/chrX_map.txt'
ref_path='./../data/chrX.fa'
reads_path='./../data/reads'
red_ranges=[149249757 , 149249868,149256127 , 149256423,149258412 , 149258580,149260048 , 149260213,149261768 , 149262007,149264290 , 149264400]
green_ranges=[149288166 , 149288277,149293258 , 149293554,149295542 , 149295710,149297178 , 149297343,149298898 , 149299137,149301420 , 149301530]


# In[16]:


sum_bit_counts=[[],[],[],[]] #ACGT
char_counts=[0,0,0,0] #ACGT
delta=50
mile=0

'''Uncomment this code to generate new Data structure files'''
# In[17]:

print('Started indexing')
with open(last_col_path,'r') as F:
    for line in F:
        for ch in line:
            if ch=='\n':
                continue
            if ch=='A':
                char_counts[0]+=1
            elif ch=='C':
                char_counts[1]+=1
            elif ch=='G':
                char_counts[2]+=1
            elif ch=='T':
                char_counts[3]+=1
            if mile%delta==0:
                for i in range(4):
                    sum_bit_counts[i].append(char_counts[i])
            mile+=1
for i in range(4):
    sum_bit_counts[i].append(char_counts[i])


# In[18]:

print('Started indexing map, needed to get reference string location from bwt column')

map_pos=[]
cnt=0
with open(chrX_map_path,'r') as C:
    for ix,line in enumerate(C):
        cnt+=len(line)
        if ix%delta==0:
            map_pos.append(cnt)
map_pos.append(cnt)
print('Finished Indexing')


print('Pickling DS')
with open('last_pre_sum.pikl','wb') as L_pkl:
    pickle.dump(sum_bit_counts,L_pkl)
with open('char_counts.pikl','wb') as cnt_pkl:
    pickle.dump(char_counts,cnt_pkl)
with open('map_pos.pikl','wb') as M_pkl:
    pickle.dump(map_pos,M_pkl)


# In[252]:


with open('last_pre_sum.pikl','rb') as L_pkl:
    sum_bit_counts=pickle.load(L_pkl)
with open('char_counts.pikl','rb') as cnt_pkl:
    char_counts=pickle.load(cnt_pkl)
with open('map_pos.pikl','rb') as M_pkl:
    map_pos=pickle.load(M_pkl)
print('Loaded pickled DS')

# In[6]:


def read_last_col(i,path,n): #reads n characters from ith character onwards from file in path ignoring newlines
    if n<=0:
        return ''
    line_no=int(i/100)
    line_pos=i-100*line_no
    f=open(path,'r')
    f.seek(line_no*101+line_pos)
    ans=''
    while (len(ans)<n):
        x=f.read(n)
        if x=='': #no more
            break
        x=x.replace('\n',"")
        ans=ans+x
        ans=ans[:n]
    f.close()
    return ans


# In[7]:


import math as m


char_to_idx={'A':0,'C':1,'G':2,'T':3}
#to begin bwt search
def init_ranges(ch):
    idx=char_to_idx[ch]
    beg=0
    for i in range(idx):
        beg+=char_counts[i]
    end=beg+char_counts[idx]-1
    return [beg,end]

def find_ranks(ch,left_col_range,path,debug=False): #return rank of first and last occurrence of ch in right column's [left_col_range[0],left_col_range[1]]
    if debug:
        print('looking for character->',ch)
    ch_idx=char_to_idx[ch]
    ans=[]
    beg,end=left_col_range
    mile_prev_beg=m.floor(beg/delta)
    mile_nxt_beg=mile_prev_beg+1
    mile_prev_end=m.floor(end/delta)
    mile_nxt_end=mile_prev_end+1
    mile_cur=mile_nxt_beg
    if debug:
        print('mile_prev',mile_prev_beg)
        print('mile_end',mile_prev_end)
    #find beginning rank
    beg_present=False
    while(mile_cur<=mile_prev_end):
        if (sum_bit_counts[ch_idx][mile_cur]!=sum_bit_counts[ch_idx][mile_nxt_beg]):
            beg_present=True
            break
        mile_cur+=1
    if debug:
        print('char present, judged from milestones',beg_present)   
    if beg_present is False:
        explicit_string=read_last_col(beg,path,min(end,mile_nxt_beg*delta)-beg+1)+read_last_col(max(beg,mile_prev_end*delta+1),path,end-max(beg,mile_prev_end*delta+1)+1)
        if ch in explicit_string:
            beg_present=True
            if debug:
                print('Char present from explicit string equality')
            
    if beg_present is False:
        return [-1,-1]
    
    beg_rank=sum_bit_counts[ch_idx][mile_prev_beg]
    rank_count_string=read_last_col(mile_prev_beg*delta+1,path,beg-1-(mile_prev_beg*delta+1)+1)
    if debug:
        print('BEG')
        print('positions rank count string',mile_prev_beg*delta+1 ,beg-1)
        print('Expected len rank count string',beg-1-(mile_prev_beg*delta+1)+1)
        print('Actual len rank count string',len(rank_count_string))
        print('Milestone count',beg_rank)
    for idx,char in enumerate(rank_count_string):
        if char==ch:
            beg_rank+=1
    if debug:
        print('Changed count',beg_rank)
            
    #edge case coinciding beg and milestone
    if beg%delta==0:
        char_at_beg=read_last_col(beg,path,1)
        if char_at_beg==ch:
            beg_rank-=1
    
    ans.append(beg_rank+1)
    
    #find last rank
    last_rank=sum_bit_counts[ch_idx][mile_nxt_end]
    rank_count_string=read_last_col(end+1,path,mile_nxt_end*delta-(end+1)+1)
    if debug:
        print('END')
        print('positions rank count string',end+1 ,mile_nxt_end*delta)
        print('Expected len rank count string',mile_nxt_end*delta-(end+1)+1)
        print('Actual len rank count string',len(rank_count_string))
        print('Milestone count',last_rank)
    for idx,char in enumerate(rank_count_string):
        if char==ch:
            last_rank-=1
    if debug:
        print('changed count',last_rank)

        
    ans.append(last_rank)
    return ans


def select(first,last,ch): #returns new range in left column from first and last ranks of character ch
    idx=char_to_idx[ch]
    beg=0
    for i in range(idx):
        beg+=char_counts[i]
    end=beg+(last-1)
    beg+=(first-1)
    return [beg,end]

def match(x,path,init_range=None,debug=False): #returns matching positions for string x, path is path of last col file
    x=x[::-1] #reverse string
    
    if init_range is None:
        left_col_range=init_ranges(x[0])
    else:
        left_col_range=init_range
    
    i=0
    not_found=False
    while(i<len(x)-1):
        if debug:
            print('current char on search string',x[i])
            print('search range',left_col_range)
            print('len range',left_col_range[1]-left_col_range[0]+1)
        nxt=i+1
        first,last=find_ranks(x[nxt],left_col_range,path,debug)
        if debug:
            print('First last for char' ,x[nxt],': ',first,', ',last)
        if first==-1:
            not_found=True
            break
        left_col_range=select(first,last,x[nxt])
        i+=1
        if debug:
            print('----------------------------------')

    if not_found:
        return []
    return left_col_range
    


# In[8]:

#reads reference location from map in O(1) time because f.seek() is O(1)
def read_map_opt(from_to,path):
    ans=[]
    f=open(path,'r')
    fr,to=from_to
    fr_beg_mile=m.floor(fr/delta)
    if fr%delta==0:
        fr_beg_mile-=1
    if fr_beg_mile==-1:
        chars_to_skip=0
    else:
        chars_to_skip=map_pos[fr_beg_mile]
    f.seek(chars_to_skip)
    delta_str=f.read(map_pos[fr_beg_mile+1]-chars_to_skip+1)
    i=fr_beg_mile*delta
    for x in delta_str:
        if i==fr-1:
            break
        chars_to_skip+=1
        if x=='\n':
            i+=1
    f.seek(chars_to_skip)
    cnt_line=0
    while(cnt_line<to-fr+1):
        L=f.readline()
        ans.append(int(L))
        cnt_line+=1
    f.close()
    return ans


# In[9]:


def generate_reverse_complement(x):
    x=x[::-1]
    reverse={'A':'T','T':'A','C':'G','G':'C'}
    ans=''
    for i,c in enumerate(x):
        ans=ans+reverse[c]
    return ans


# In[10]:

#replaces N to A
def clean_read(x):
    ans=''
    for c in x:
        if c!='A' and c!='T' and c!='G' and c!='C':
            c='A'
        ans=ans+(c)
    return ans

#wrapper function to match x with references and return all occurrences
def match_read(x): 
    x=clean_read(x)
    ans_range=match(x,last_col_path)
    if len(ans_range)==0:
        return [-1]
    ref_position=read_map_opt(ans_range,chrX_map_path)
    return ref_position


# In[11]:

#
def read_reference(i,path,n):
    if i<0:
        return ''
    f=open(path,'r')
    base=6
    if n<=0:
        return ''
    line_no=int(i/100)
    line_pos=i-100*line_no
    f.seek(base+line_no*101+line_pos)
    ans=''
    while (len(ans)<n):
        x=f.read(n)
        if x=='': #no more
            break
        x=x.replace('\n',"")
        ans=ans+x
        ans=ans[:n]
    f.close()
    return ans


# In[12]:

#match with relaxed error of 2 mismatches
def matches(x,y):
    mismatch=0
    if(len(x)==len(y)):
        return False
    for i,xx in enumerate(x):
        if y[i]!=xx:
            mismatch+=1
    return True if mismatch<=2 else False

#breaks the string into 3 parts and 1 of them should match. This is taken as the reference for others.
def mismatched_match(read):
    L=len(read)
    read_a=read[:int(L/3)]
    read_b=read[int(L/3):2*int(L/3)]
    read_c=read[2*int(L/3):]
#     print(len(read_a),len(read_b),len(read_c))
#     print(read_a)
#     print(read_b)
#     print(read_c)
    assert(len(read_a)+len(read_b)+len(read_c)==len(read))
    match_a=match_read(read_a)
    match_b=match_read(read_b)
    match_c=match_read(read_c)
    ans=[]
#     print(match_a,match_b,match_c)
    for pos_a in match_a:
        if pos_a==-1:
            continue
        str_to_match=read_reference(pos_a,ref_path,L)
#         print(str_to_match,read)
        if matches(str_to_match,read):
            ans.append(pos_a)
            
    for pos_b in match_b:
        if pos_b==-1:
            continue
        str_to_match=read_reference(pos_b-len(read_a),ref_path,L)
        if matches(str_to_match,read):
            ans.append(pos_b-len(read_a))
            
    for pos_c in match_c:
        if pos_c==-1:
            continue
        str_to_match=read_reference(pos_c-len(read_a)-len(read_b),ref_path,L)
        if matches(str_to_match,read):
            ans.append(pos_c-len(read_a)-len(read_b))
    if len(ans)==0:
        return [-1]
    return ans   


# In[13]:

#get exon from matched positions for 1 read
def get_exon(matched_pos,r,g):
    r_matched=[] #r exons matched for this read
    g_matched=[] #g exons matched for this read
    for z in matched_pos:
        i=0
        while( i < len(red_ranges)):
            nxt=i+1
            if red_ranges[i]<=z and red_ranges[nxt]>=z:
                r_matched.append(int(i/2))
            i+=2
            
        while( i < len(green_ranges)):
            nxt=i+1
            if green_ranges[i]<=z and green_ranges[i]>=z:
                g_matched.append(int(i/2))
            i+=2
    k=len(r_matched)
    for x in r_matched:
        r[x]+=(1/k) #assigning equal weights
    if k!=0:
        return (r,g) #because assign all green matched positions to red matched by convention
    k=len(g_matched)
    for x in g_matched:
        g[x]+=(1/k)
    return (r,g)


# In[35]:


import time
c=0
r=[0,0,0,0,0,0] #count of matches red and green exons 
g=[0,0,0,0,0,0]
unmatched=0
with open(reads_path,'r') as F:
    start=time.time()
    for ix,line in enumerate(F):
        #uncomment these 2 lines to run on whole read set.
        if ix<2906996 or ix>=2941395:
            continue
        #remove new line
        read=line[:-1]
        read=clean_read(read)
        #check for direct match
        match_pos=match_read(read)
        if match_pos[0]==-1:
            #check for reverse complement match
            match_pos=match_read(generate_reverse_complement(read))
        if match_pos[0]==-1:
            match_pos=mismatched_match(read)
            if match_pos[0]==-1:
                match_pos=mismatched_match(generate_reverse_complement(read))
        r,g=get_exon(match_pos,r,g)
        if match_pos[0]==-1:
            unmatched+=1
        c+=1
        if c%100==0:
            end_c=time.time()
            print(ix,end_c-start)
            start=end_c
            # print('Red green counts',r,g)

from scipy.stats import binom
#returns probability of observing counts given the configuration
def prob(config,counts):
    print(counts)
    ans=1.0
    for ix in range(len(counts[0])):
        if ix==0 or ix==5:
            continue
        n_r=int(counts[0][ix])
        n_g=int(counts[1][ix])
        p=(config[ix]/100)/((config[ix]/100)+1)
        ans=ans* binom.pmf(n_r,n_r+n_g,p)
    return ans


# In[47]:


configs=[[0,50,50,50,50,0],[0,100,100,0,0,0],[0,33,33,100,100,0],[0,33,33,33,100,0]]
p=[]
for config in configs:
    p.append(prob(config,[r,g]))


# In[48]:

print('Probabilities for each distribution ,0 to 4 as in slide, is,',p)
# print(match_read('GAGGACAGCACCCAGTCCAGCATCTTCACCTACACCAACAGCAACTCCACCAGAGGTGAGCCAGCAGGCCCGTGGAGGCTGGGTGGCTGCACTGGGGGCCA'))