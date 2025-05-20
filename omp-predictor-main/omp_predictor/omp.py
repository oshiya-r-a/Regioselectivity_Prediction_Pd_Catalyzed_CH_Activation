def main():
    import argparse
    import os

    parser = argparse.ArgumentParser(description="Predict the regioselectivity from a .mol structure using a pre-trained ML models.")
    parser.add_argument("input_folder", type=str, help="Folder containing .mol files")
    parser.add_argument("--lf", action="store_true", help="Use largefragment encoder (default is standard)")
    parser.add_argument("--logreg", action="store_true", help="Use Logistic Regression model (default is SVM)")
    args = parser.parse_args()

    from . import fingerprint, encoder

    fingerprint_type = "lf" if args.lf else "standard"
    predictor = "logreg" if args.logreg else "svm"

    from rdkit import Chem
    import joblib

    mol_files = [os.path.join(args.input_folder, f) for f in os.listdir(args.input_folder) if f.endswith(".mol")]
    mols = [Chem.MolFromMolFile(file_name) for file_name in mol_files]

    morgan_fingerprint = fingerprint.MorganFingerprint(
        nbits=1024,
        radius=2,
        features=True,
        chirality=True
    )

    if fingerprint_type == "lf":
        smaller_fragments_morgan_fingerprint = fingerprint.MorganFingerprint(
            nbits=1024,
            radius=1,
            features=True,
            chirality=True
        )
        bit_vectors = encoder.largefragment_encoder(mols, morgan_fingerprint, smaller_fragments_morgan_fingerprint)
        model_file = "lf_svm_model.pkl" if predictor == "svm" else "lf_lg_model.pkl"

    else:
        bit_vectors = encoder.standard_encoder(mols, morgan_fingerprint)
        model_file = "svm_model.pkl" if predictor == "svm" else "lg_model.pkl"

    model_path = os.path.join(os.path.dirname(__file__), model_file)
    model = joblib.load(model_path)

    class_names = ["meta", "non", "ortho", "para"]

    import numpy as np
    predictions = model.predict(bit_vectors)
    probabilities = np.max(model.predict_proba(bit_vectors), axis=1)
    for file_name, prediction, probability in zip(mol_files, predictions, probabilities):
        if prediction != 1:
            if probability > 0.6:
                print(f" {os.path.basename(file_name)} is a {class_names[prediction]}-director with {int(probability*100)}% certainty.")
            else:
                print(f" {os.path.basename(file_name)} is a non-director.")
        
        else:
            print(f" {os.path.basename(file_name)} is a {class_names[prediction]}-director.")


if __name__ == "__main__":
    main()