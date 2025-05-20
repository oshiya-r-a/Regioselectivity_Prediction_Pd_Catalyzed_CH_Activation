from rdkit.Chem import AllChem
import numpy as np

class Fingerprint:
    
    # This class is is equvalent to a data class which contains information about morgan fingerprint
    
    fingerprints = []
    
    def __init__(self, mol, fingerprint, on_bits, bit_info):
        self.mol = mol # molecule from which the fingerprint is generated
        self.fingerprint = fingerprint # bit vector
        self.on_bits = on_bits # set of indices in bit vector having a value of 1
        self.bit_info = bit_info # fragments encoded in the bits having value of 1
        Fingerprint.fingerprints.append(self) # Keeps a record of each fingerprint created


class MorganFingerprint:
    
    # This class contains the parameters for Morgan fingerprints and featurize funtion
    
    mols = []
    fps = []
    
    def __init__(self, nbits:int, radius:int, features:bool, chirality:bool):
        self.nbits = nbits
        self.radius = radius
        self.features = features
        self.chirality = chirality
        
    def featurize(self, mol):
        MorganFingerprint.mols.append(mol)
        bit_info = {}
        fp = AllChem.GetMorganFingerprintAsBitVect(
                mol,
                self.radius,
                nBits = self.nbits,
                useFeatures = self.features,
                useChirality = self.chirality,
                bitInfo = bit_info
            )
        
        MorganFingerprint.fps.append(fp)
        
        return Fingerprint(mol=mol,
                           fingerprint=np.array(fp),
                           on_bits=sorted(set(fp.GetOnBits())),
                           bit_info=bit_info
                          )