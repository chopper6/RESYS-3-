from math import log, pow, factorial as fct, floor
import random as rd, sys
from scipy.special import comb

# Pr should be a dictionary of probabilities, as built in pr.py

debug = False
ROUND=12
sys.setrecursionlimit(10000)

########################################### WHOLE #################################################

def H(Pr, logbase=2):
    num_inst = len(Pr)
    H=0
    for event in Pr:
    	p = Pr[event]
    	#print(event,p,-1*p*log(p,2))
    	if p!=0: H+= -1*p*log(p,logbase)

    # use logbase = #bits???
    return H


def Ee_empirical(n=3940, mMax=40, reps=100):
	Ee = [0 for i in range(mMax)]
	for r in range(reps):
		H=[0 for i in range(mMax)]
		for m in range(mMax):
			bins = {}
			for i in range(n):
				event = int(rd.random()*pow(2,m))
				if event in bins.keys():
					bins[event] += 1
				else:
					bins[event] = 1
			summ = sum(bins[b] for b in bins.keys())
			for b in bins.keys():
				pr=bins[b]/summ
				if pr!=0:
					H[m] -= pr * log(pr,2)
			#print(m,[bins[b]/summ for b in bins.keys()],H[m],'\n')

		for i in range(mMax):
			Ee[i] += H[i]
		curr_Ee = [Ee[i]/(r+1) for i in range(mMax)]
		#if r%10==0: print('Ee at iter',r,':',curr_Ee)

	Ee = [Ee[i]/reps for i in range(mMax)]
	print("\n\nEe Final = ",Ee)
	return Ee

def expected_entropy(m,n):
	lng, summ = int(pow(2,n)),m # 
	Ee = gen_strings(lng,summ,summ,[])
	return Ee


def expected_entropy_memOF(m,n):
	lng, summ = int(pow(2,n)),m # 

	S = gen_strings(lng,summ,summ)
	print("Done generating string")
	Ee = 0
	for s in S:
		Ee += combos_entropy(s)
	return Ee

def combos_entropy(string):
	# where string is an integer string corresp to a binary matrix
	# calc #combos = #int_perms * #matrix_perms
	int_perms = round(fct(sum(string)), ROUND)
	for s in string:
		int_perms = int_perms // round(fct(s),ROUND)
	unq_s={}
	for s in string:
		if s in unq_s.keys():
			unq_s[s] += 1
		else: 
			unq_s[s] = 1
	matrix_perms = round(fct(len(string)),ROUND)#(sum(unq_s))
	H=0
	for s in unq_s.keys():
		matrix_perms = matrix_perms // round(fct(unq_s[s]),ROUND)

	for s in string:
		pr = s/sum(string)
		assert(pr>=0 and pr<=1)
		if pr!=0: H -= pr*log(pr,2)

	total2 = pow(2,sum(string)*log(len(string),2))
	#print(string,int_perms,matrix_perms,total2,H)
	prob_perm = int_perms*matrix_perms/total2
	assert(prob_perm>=0 and prob_perm<=1)

	return prob_perm*H


def gen_strings(lng,summ, prev, partstring):
	# init: lng = 2^m, summ = n
	# note that first string
	if lng ==0 and summ>0: return 0
	if lng==0: return combos_entropy(partstring)
	if summ==0: return combos_entropy(partstring+[0 for i in range(lng)])
	T = 0
	for i in range(1,min(prev,summ)+1):
		T += gen_strings(lng-1,summ-i,i,partstring+[i])
	#print('gen_str(',lng,summ,') returns:',S)
	return T





################################## UNUSED FNS BELOW ######################################


def H_cond(Pr,a,b):
	# H(a|b)
    num_inst = len(Pr)
    H = sum(h_cond(Pr[i],a,b) for i in range(num_inst)) / num_inst
    return H

def Info(Pr,a,b):
    num_inst = len(Pr)
    I = sum(info(Pr[i],a,b)for i in range(num_inst)) / num_inst
    return I




########################################### POINTWISE ############################################

log_base = 2

def h(pr,a,logbase=2):
    # h(a)
    return -1*log(pr[a], logbase)

def h_cond(pr,a,b):
    # h(a|b)
    return -1*log(pr[str(a) + str(b)]/pr[b], log_base)

def info(pr,a,b):
    # i(a,b)
    return log(pr[str(a) + str(b)]/(pr[a]*pr[b]), log_base)


if __name__ == "__main__":
	m, n = 20,10 #3934,10
	print(expected_entropy(m,n))

	if False:
		Ee = []
		for i in range(n+1):
			Ee += [expected_entropy(m,i)]
			print(i,Ee)

		print("Final expected entropy array:\n",Ee)

