def eval_jaccard(predFile,groundTruthFile):

	pred = []
	with open(predFile) as f:
		lines = f.read().splitlines()
		for line in lines:
			split = line.split()
			community = []
			for x in split[1:len(split)]:
				community.append(int(x))
			pred.append(community)

	
	actual = []

	with open(groundTruthFile) as f:
		lines = f.read().splitlines()
		for line in lines:
			split = line.split()
			community = []
			for x in split[1:len(split)]:
				community.append(int(x))
			actual.append(community)


	score1 = 0.0
	score2 = 0.0


	for i in pred:
		max_i = 0
		for j in actual:
			jaccard_sim = len(set(i).intersection(set(j)))*1.0/ len(set(i).union(set(j)))
			if max_i < jaccard_sim:
				max_i = jaccard_sim
		score1 += max_i


	score1 /= (2.0 * len(pred))

	for i in actual:
		max_i = 0
		for j in pred:
			jaccard_sim = len(set(i).intersection(set(j)))*1.0/ len(set(i).union(set(j)))
			if max_i < jaccard_sim:
				max_i = jaccard_sim
		score2 += max_i


	score2 /= (2.0 * len(actual))

	score = score1 + score2
	return score


def eval_f1(predFile, groundTruthFile):
	

	pred = []
	with open(predFile) as f:
		lines = f.read().splitlines()
		for line in lines:
			split = line.split()
			community = []
			for x in split[1:len(split)]:
				community.append(int(x))
			pred.append(community)

	
	actual = []

	with open(groundTruthFile) as f:
		lines = f.read().splitlines()
		for line in lines:
			split = line.split()
			community = []
			for x in split[1:len(split)]:
				community.append(int(x))
			actual.append(community)


	score1 = 0.0
	score2 = 0.0

	for i in pred:
		max_i = 0
		for j in actual:
			p =len(set(i).intersection(set(j)))*1.0/ len(set(i))
			r =len(set(i).intersection(set(j)))*1.0/ len(set(j))
			if p > 0.0 and r >0.0:
				f1 = 2*p*r/(p+r)
			else:
				f1 = 0.0
			if max_i < f1:
				max_i = f1
		score1 += max_i


	score1 /= (2.0 * len(pred))

	for i in actual:
		max_i = 0
		for j in pred:
			p =len(set(i).intersection(set(j)))*1.0/ len(set(i))
	        r =len(set(i).intersection(set(j)))*1.0/ len(set(j))
	        if p>0.0 and r > 0.0:
	        	f1 = 2*p*r/(p+r)
	        else:
	        	f1 = 0.0
	        if max_i < f1:
	        	max_i = f1
		score2 += max_i


	score2 /= (2.0 * len(actual))

	score = score1 + score2
	return score

avg_jac = 0.0
avg_f1 = 0.0
print "jaccard", "f1"
with open('facebook/filenames.txt') as f: 
	lines = f.read().splitlines()
	for f in lines:
		p1 = "facebook/"+f
		p2 = "pred_net/"+f
		# p2 = "Girvan_Results/"+f
		print eval_jaccard(p1,p2), eval_f1(p1,p2)
		avg_f1 += eval_f1(p1,p2)
		avg_jac += eval_jaccard(p1,p2)

print "Avg_jac", "Avg_f1"		
print avg_jac/10, avg_f1/10