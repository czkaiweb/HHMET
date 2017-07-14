import sys

dataname = sys.argv[1]

MF_l=[ i*50+200 for i in range(10) ]
MY_l=[ i*100+600 for i in range(8)]
gYqq_l=[ 0.05+i*0.05 for i in range(10)]
gYgg_l=[ 0.01+i*0.01 for i in range(20)]

if 'qq' in dataname:
#	fout=open("run_dijet_qq.sh","w+")
	g_coupling = gYqq_l
if 'gg' in dataname:
#	fout=open("run_dijet_gg.sh","w+")
	g_coupling= gYgg_l

foutname = "run_script/run_"+dataname

fout=open("analysis.sh","w+")
fout.write("#! /bin/bash\n")

for i in range(len(MY_l)):
	for j in range(len(MF_l)):
		for k in range(len(g_coupling))	:
                        if 'qq' in dataname:
#				filename=foutname+"gYqq_"+str(g_coupling[k])+"_MY_"+str(int(MY_l[i]))+"_MF_"+str(int(MF_l[j]))+".sh"
				cardname='gYqq'+str(g_coupling[k])+'MY'+str(int(MY_l[i]))+'MF'+str(int(MF_l[j])) #Full
			if 'gg' in dataname:
#				filename=foutname+"_gYgg_"+str(int(100*g_coupling[k]))+"_MY_"+str(int(MY_l[i]))+"_MF_"+str(int(MF_l[j]))+".sh"
				cardname='gYgg'+str(g_coupling[k])+'MY'+str(int(MY_l[i]))+'MF'+str(int(MF_l[j])) #Full
			fout.write("./run_analysis.sh %s\n"%(cardname))

fout.close()
