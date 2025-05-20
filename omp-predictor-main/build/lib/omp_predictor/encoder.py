def standard_encoder(mols, fingerprint):
    X = []
    for mol in mols:
        X.append(fingerprint.featurize(mol).fingerprint.tolist())
    return X

def largefragment_encoder(mols, fingerprint_large_radius, fingerprint_small_radius):
    X = []
    for mol in mols:
        fp_all_fragments = fingerprint_large_radius.featurize(mol).fingerprint.tolist()
        fp_smaller_fragments = fingerprint_small_radius.featurize(mol).fingerprint.tolist()
        
        smaller_fragment_indexes = []
        for i in range(len(fp_smaller_fragments)):
            if fp_smaller_fragments[i] == 1:
                smaller_fragment_indexes.append(i)
            i += 1
        
        for index in smaller_fragment_indexes:
            fp_all_fragments[index] = 0
            
        X.append(fp_all_fragments)
        
    return X