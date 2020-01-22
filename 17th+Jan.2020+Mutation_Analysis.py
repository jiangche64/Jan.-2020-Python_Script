
# coding: utf-8

# # GOAL: Mutation Files

# ## STEP 1 Exon Extraction

# In[1]:

Exon_list= []
with open ('gencode.vM21.annotation.gtf') as file:
    for line in file:
        if line.startswith('#'):
            continue
        debris = line[:-1].split('\t')
        cata = debris[2]
        if cata =='exon':
            chro = debris[0]
            start= int(debris[3])
            stop = int(debris[4])
            ge_id= debris[8].split(';')[0][9:-1]
            Exon_list.append([chro, ge_id, start, stop])


# In[4]:

Exon_list #BROVA


# ## STEP 2 Mutation Info Extraction

# In[2]:

B1 =[]
with open ('B1.vcf') as file:
    for line in file:
        if line.startswith('#'):
            continue
        debris = line.split('\t')
        qual = float(debris[5])
        if qual >= 99:
            chro = debris[0]
            posi = int(debris[1])
            ref = debris[3]
            mut = debris[4]
            B1.append([chro, posi, ref, mut])


# In[9]:

B1


# ## STEP 3 Mutation Analysis

# In[3]:

B1_mu = {}
atcg = ['A','T','C','G']
for x in atcg:
    B1_mu[x] ={}
    for y in atcg:
        B1_mu[x][y] = 0

for x in atcg:
    for y in atcg:
        for b in range(0, len(B1)):
            if B1[b][2] == x and B1[b][3] == y:
                B1_mu[x][y] += 1


# In[4]:

B1_mu


# ## Multi-files in defination

# In[5]:


def Bfile_mu_1(infile, B_mu_dict): # infile = str , Blist = []
    Blist = []
    with open (infile) as file:
        for line in file:
            if line.startswith('#'):
                continue
            debris = line.split('\t')
            qual = float(debris[5])
            if qual >= 99:
                chro = debris[0]
                posi = int(debris[1])
                ref = debris[3]
                mut = debris[4]
                Blist.append([chro, posi, ref, mut])

    atcg = ['A','T','C','G']
    for x in atcg:
        B_mu_dict[x] ={}
        for y in atcg:
            B_mu_dict[x][y] = 0

    for x in atcg:
        for y in atcg:
            for b in range(0, len(Blist)):
                if Blist[b][2] == x and Blist[b][3] == y:
                    B_mu_dict[x][y] += 1


# In[6]:

B1_mu_dict = {}
Bfile_mu_1('B1.vcf', B1_mu_dict)


# In[7]:

B1_mu_dict # Congratulations!


# In[8]:

B2_mu_dict = {}
Bfile_mu_1('B2.vcf', B2_mu_dict)


# In[9]:

B2_mu_dict


# # Splendid
