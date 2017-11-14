def consensus(sequences, qualities):
    consensus_quality = ''
    consensus = ''
    for i in range(max((len(s) for s in sequences))):
        votes = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for idx, seq in enumerate(sequences):
            if len(seq) <= i:
                continue
            votes[seq[i]] += math.log10(1-phred(qualities[idx][i]))
            for b in votes:
                if b == seq[i]:
                    continue
                votes[b] += math.log10(phred(qualities[idx][i]))
        for v in votes:
            votes[v] = math.pow(10, votes[v])
        norm_fact = sum(votes.values())
        for v in votes:
            votes[v] = votes[v]/norm_fact
        #print(votes)
        consensus += max(votes.keys(), key=(lambda k: votes[k]))
        consensus_quality += qual(1-votes[consensus[-1]])
    return consensus, consensus_quality
