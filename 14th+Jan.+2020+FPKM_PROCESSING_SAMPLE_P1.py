
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

# # Double files check

# In[ ]:

#DATA from 'wll_mito_gene_matrix.txt'

symble	UT_1_SL_P1	UT_2_SL_P1	UT_3_SL_P1	DMM_20_1_SL_P1	DMM_20_2_SL_P1	DMM_20_3_SL_P1	3NPA_100_1_SL_P1	3NPA_100_2_SL_P1	3NPA_100_3_SL_P1	Succinate10_1_SL_P1	Succinate10_2_SL_P1	Succinate10_3_SL_P1	SL_2i_1_P1	SL_2i_2_P1	SL_2i_3_P1	UT_1_SL_P3	UT_2_SL_P3	UT_3_SL_P3	DMM_20_1_SL_P3	DMM_20_2_SL_P3	DMM_20_3_SL_P3	3NPA_100_1_SL_P3	3NPA_100_2_SL_P3	3NPA_100_3_SL_P3	Succinate10_1_SL_P3	Succinate10_2_SL_P3	Succinate10_3_SL_P3	SL_2i_1_P3	SL_2i_2_P3	SL_2i_3_P3
Gnai3	57.925	63.149	57.369	62.259	70.327	50.658	54.859	61.84	55.274	60.865	64.615	49.822	48.416	47.36	43.798	59.935	50.274	50.321	50.063	35.477	51.479	50.461	42.083	51.505	55.426	50.783	53.762	38.908	28.402	41.571
Pbsn	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
Cdc45	7.671	8.474	7.926	6.702	7.727	7.949	7.1	8.082	7.966	7.149	9.126	7.518	7.358	9.261	8.265	8.839	7.927	8.192	10.242	11.227	9.16	8.734	8.865	8.131	10.419	8.958	7.809	12.558	11.745	9.83


# In[58]:

Gnai3 = ''
for x in [0,5,10, 1,6,11]:
    Gnai3 += str(symble_FPKM['Gnai3'][x])+'\t'
Cdc45 = ''
for x in [0,5,10, 1,6,11]:
    Cdc45 += str(symble_FPKM['Cdc45'][x])+'\t'


# In[59]:

Gnai3 #DATA from JC_FILE


# In[60]:

Cdc45 #DATA from JC_FILE


# In[ ]:

# Serum + Lif 	 E14 	 P1 
symble	rep1_UT	rep1_Treat_D	rep1_Treat_N	rep1_Treat_S	rep1_Treat_i	rep2_UT	rep2_Treat_D	rep2_Treat_N	rep2_Treat_S	rep2_Treat_i	rep3_UT	rep3_Treat_D	rep3_Treat_N	rep3_Treat_S	rep3_Treat_i
Gnai3	57.76488596658472	62.08704504468488	54.70802376279869	60.69732024984706	48.28209737815031	62.97496484988169	70.13305318269515	61.66942263649407	64.43644038250729	47.22953184114676	57.2108731693898	50.517967662344105	55.12129507080678	49.684729982290776	43.677476243025346	
Cdc45	7.616531918041945	6.654462922127366	7.049669970276409	7.098313245765594	7.306208109401504	8.414280592295937	7.67265576982131	8.024909710868853	9.061471199484414	9.195107356348695	7.869873323719514	7.89252661845749	7.909408272128963	7.464215543111708	8.206396251814187	
H19	0.6542296410028855	0.984520290534543	0.47594714934358184	0.671763943966716	1.1223572059309785	0.48747303345641196	0.7668408242696869	0.26395301915362984	0.48924342328753184	1.0323882500704453	0.5899137463780332	1.4690716849014	0.7251931734811572	1.3456479563258825	1.589423555385416	

