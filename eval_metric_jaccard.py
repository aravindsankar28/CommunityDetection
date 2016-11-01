import sys

predFile = sys.argv[1]

pred = []
with open(predFile) as f:
	lines = f.read().splitlines()
	for line in lines:
		split = line.split()
		community = []
		for x in split[1:len(split)]:
			community.append(int(x))
		pred.append(community)

groundTruthFile = sys.argv[2]
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

print "Score = ",score
