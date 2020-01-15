
# coding: utf-8

# # FPKM (a_whole_processing_file)
# Step 1 Assemble all the files (B1 to B30)
# Step 2 Remove the nonsence value( [0]*15 )
# Step 3 Caculate the len(exons) (from 'gencode.vM21.annotation.gtf')
# Step 4 Caculate the Total_Reads_per_sample
# Step 5 Result = Reads/ len(Exon)/ Total_Reads_per_sample*1e9

# # Step1
# Target: I want to pack all those data into a dictionary
# Have a test on silgle/double file(s), them in use into a function

# In[5]:

out_dic = {}
with open ('B (1).txt') as file:
    for line in file:
        if line.startswith('_'):
            continue
        debris = line[:-1].split('\t')
        name = debris[0]
        expr = int(debris[1])
        out_dic[name]= [expr]

with open ('B (2).txt') as file:
    for line in file:
        if line.startswith('_'):
            continue
        debris = line[:-1].split('\t')
        name = debris[0]
        expr = int(debris[1])
        out_dic[name].extend([expr])
        


# In[4]:

out_dic


# In[6]:

def dict_processing1(infile , out_dict): #infile = str , out_dic = dict
#     out_dic = {}
    with open (infile) as file:
        for line in file:
            if line.startswith('_'):
                continue
            debris = line[:-1].split('\t')
            name = debris[0]
            expr = int(debris[1])
            out_dict[name]= [expr]

def dict_processing2(infile , out_dict):  #infile = str , out_dic = dict(from dict_processing1)
    with open (infile) as file:
        for line in file:
            if line.startswith('_'):
                continue
            debris = line[:-1].split('\t')
            name = debris[0]
            expr = int(debris[1])
            out_dict[name].extend([expr])
        


# In[7]:

expr_assemblage = {}

dict_processing1('B (1).txt', expr_assemblage)

for x in (list(range(2,6)) + list(range(11,16)) + list(range(21,26))):
    file_name = 'B ('+ str(x) + ').txt'
    dict_processing2(file_name, expr_assemblage)
    
# dict_processing2('B (2).txt', expr_assemblage)
# dict_processing2('B (3).txt', expr_assemblage)
# dict_processing2('B (4).txt', expr_assemblage)
# dict_processing2('B (5).txt', expr_assemblage)
# dict_processing2('B (6).txt', expr_assemblage)


# # Step2
# Target: Remove the items with [0]*15
# Attention: Set up a new can for storage

# In[8]:

new_assemblage = {}
for x in expr_assemblage:
    if expr_assemblage[x].count(0) != 15: #that means an item in a list, not an item in a number! See the example below
        new_assemblage[x] = expr_assemblage[x]


# In[53]:

new_assemblage


# # Step3
# Target: 
# exons in a gene

# In[9]:

exon_len = {}
with open('../Python_2019_Dec/mut/gencode.vM21.annotation.gtf') as file:
    for line in file:
        if line.startswith('#'):
            continue
        debris = line[:-1].split('\t')
        types = debris[2]
        if types == 'exon':
            ge_id = debris[8].split(';')[0][9:-1] # gene_id in use, which is cordinate with those formal files
            start = int(debris[3])
            stop  = int(debris[4])
            length= abs(start-stop) +1 # Do take care
            if ge_id in exon_len:
                exon_len[ge_id] += length # multi-exons may found in a single gene, so combine them together
            else:
                exon_len[ge_id] = length


# In[46]:

exon_len


# # Step4
# Target: Do a summation of all the reads per column
# Tips: Function 'range' been use for index & name

# In[10]:

sum_reads = {}
for x in range(15):
    sum_reads[x] = 0
    for y in new_assemblage:
        sum_reads[x] += new_assemblage[y][x] # A bit calculating redundant! Do need a simplification.


# In[55]:

sum_reads


# # Step5
# Computing the FPKM

# In[11]:

value_FPKM = {}
for x in new_assemblage:
    value_FPKM[x] = []
    for a in range(15):
        value_FPKM[x].append(new_assemblage[x][a]/exon_len[x]/sum_reads[a]*1e9)


# In[ ]:

value_FPKM


# # Step6
# Convert those gene_ids to symbles
# And Export the results

# In[15]:

symble_dict = {}
with open ('mm10_gene_exon_long.txt') as file:
    for line in file:
        if not line.startswith('gene'):
            piece = line[:-1].split('\t')
            ge_id = piece[0]
            symble= piece[3]
            symble_dict[ge_id] = symble

# for ge_id , symble in symble_dict.items():
#     symble_FPKM[symble] = value_FPKM[ge_id]

symble_FPKM = {}
for ge_id, value in value_FPKM.items():
    symble_FPKM[symble_dict[ge_id]] = value_FPKM[ge_id]


# In[14]:

symble_dict


# In[16]:

symble_FPKM


# In[42]:

treatments = ['UT', 'Treat_D', 'Treat_N', 'Treat_S', 'Treat_i']
with open('14th Jan. 2020 FPKM_PROCESSING_SAMPLE_P1.txt','w') as file:
    file.write('# Serum + Lif \t E14 \t P1 \n''symble')
    for x in range(1, 4):
        for t in treatments:
            file.write('\t'+'rep' +str(x)+'_'+ str(t))
                   
                   
    for x in symble_FPKM:
        file.write('\n' + str(x) +'\t' )
        for y in range(15):
            file.write(str(symble_FPKM[x][y]) + '\t')


# # BRAVOO!
