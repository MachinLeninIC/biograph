from Bio import pairwise2
import operator

def get_alignment(sequenceA, sequenceB, featuresA, featuresB):
    """
    Yields a generator of features for aligned sequences.

    Parameters
    ----------
    sequenceA: string, biopython sequence objects or list.
        The first sequence.
    sequenceB: string, biopython sequence objects or list.
        The second sequence.
    featuresA: iterable
        Features for the first sequence. Must have the same length as sequenceA.
    featuresB: iterable
        Features for the second sequence. Must have the same length as sequenceB.

    Returns
    -------
    featureA: element of featuresA aligned with element of featuresB
    featureB: element of featuresB aligned with element of featuresA 
    """
    alignment = max(pairwise2.align.globalxx(sequenceA, sequenceB), key=operator.itemgetter(1))
    alignedA, alignedB, _, _, _ = alignment

    featuresA = iter(featuresA)
    featuresB = iter(featuresB)

    for x,y in zip(alignedA, alignedB):
        Fa = next(featuresA) if x != '-' else None
        Fb = next(featuresB) if y != '-' else None

        if Fa and Fb:
            yield Fa, Fb

def join_conservation_data(sequence, features_dict, conservation_file):
    """
    Joins feature_dict with the conservation features based on the sequence
    alignment. Features loaded are score, color, score confidence interval,
    color confidence interval and residue variety.

    Parameters
    ----------
    sequence: string, biopython sequence objects or list.
        The sequence to be used in alignment.
    features_dict: dictionary
        Dictionary mapping residues to features.
    conservation_file: string
        Full path to consurf file to be used for features.

    Returns
    -------
    features_dict: dictionary
        Same dictionary as before, but with the new features if the residue
        had a match in the consurf file.
    """
    with open(conservation_file, 'r') as ifile:
        lines = [x.rstrip('\n') for x in ifile.readlines()]
        lines = [x.split('\t') for x in lines]
        lines = [x for x in lines if len(x) == 14]

    features = features_dict.keys()
    cons_sequence = ''.join([line[1].lstrip() for line in lines])
    for res_id, line in get_alignment(sequence, cons_sequence, features, lines):
        features_dict[res_id]["score"] = float(line[3])
        features_dict[res_id]["color"] = line[5]
        features_dict[res_id]["score_confidence_interval"] = line[6]
        features_dict[res_id]["color_confidence_interval"] = line[9]
        features_dict[res_id]["residue_variety"] = line[-1]
    return features_dict
