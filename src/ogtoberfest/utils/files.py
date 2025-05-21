import os
from ogtobergest.utils import util
    
    
def process_fasta_files(fasta_dir, species_ids_fn, sequence_ids_fn):
    fasta_extensions = {"fa", "faa", "fasta", "fas", "pep", "fna"}

    if not os.path.exists(fasta_dir):
        print(f"\nDirectory does not exist: {fasta_dir}")
        util.Fail()

    files_in_directory = sorted([
        f for f in os.listdir(fasta_dir) if os.path.isfile(os.path.join(fasta_dir, f))
    ])
    fasta_files = []
    excluded_files = []

    for f in files_in_directory:
        if len(f.rsplit(".", 1)) == 2 and f.rsplit(".", 1)[1].lower() in fasta_extensions and not f.startswith("._"):
            fasta_files.append(f)
        else:
            excluded_files.append(f)

    if len(excluded_files) != 0:
        print("\nWARNING: Files have been ignored as they don't appear to be FASTA files:")
        for f in excluded_files:
            print(f)
        print("OrthoFinder expects FASTA files to have one of the following extensions: %s" % (", ".join(fastaExtensions)))
    
    if len(fasta_files) < 2:
        print("ERROR: At least two species are required")
        if len(fasta_files) == 0:
            print(f"\nNo fasta files found in supplied directory: {fasta_dir}")
        util.Fail()

    iSeq = 0
    iSpecies = 0

    newSpeciesIDs = []

    with open(sequence_ids_fn, 'a') as sequence_ids_file, \
    open(species_ids_fn, 'a') as species_ids_file:
        for fastaFilename in originalFastaFilenames:
            newSpeciesIDs.append(iSpecies)
            outputFasta = open(files.FileHandler.GetSpeciesFastaFN(iSpecies, qForCreation=True), 'w')
            fastaFilename = fastaFilename.rstrip()
            species_ids_file.write("%d: %s\n" % (iSpecies, fastaFilename))
            baseFilename, extension = os.path.splitext(fastaFilename)
            mLinesToCheck = 100
            qHasAA = False
            with open(fastaDir + os.sep + fastaFilename, 'r') as fastaFile:
                for iLine, line in enumerate(fastaFile):
                    if line.isspace(): continue
                    if len(line) > 0 and line[0] == ">":
                        newID = "%d_%d" % (iSpecies, iSeq)
                        acc = line[1:].rstrip()
                        if len(acc) == 0:
                            print("ERROR: %s contains a blank accession line on line %d" % (fastaDir + os.sep + fastaFilename, iLine+1))
                            util.Fail()
                        sequence_ids_file.write("%s: %s\n" % (newID, acc))
                        outputFasta.write(">%s\n" % newID)    
                        iSeq += 1
                    else:
                        line = line.upper()    # allow lowercase letters in sequences
                        if not qHasAA and (iLine < mLinesToCheck):
#                            qHasAA = qHasAA or any([c in line for c in ['D','E','F','H','I','K','L','M','N','P','Q','R','S','V','W','Y']])
                            qHasAA = qHasAA or any([c in line for c in ['E','F','I','L','P','Q']]) # AAs minus nucleotide ambiguity codes
                        outputFasta.write(line)
                outputFasta.write("\n")
            if (not qHasAA) and (not q_dna):
                qOk = False
                print("ERROR: %s appears to contain nucleotide sequences instead of amino acid sequences. Use '-d' option" % fastaFilename)
            iSpecies += 1
            iSeq = 0
            outputFasta.close()
        if not qOk:
            util.Fail()

    if len(originalFastaFilenames) > 0: outputFasta.close()
    speciesInfoObj.speciesToUse = speciesInfoObj.speciesToUse + newSpeciesIDs
    speciesInfoObj.nSpAll = max(speciesInfoObj.speciesToUse) + 1      # will be one of the new species
    
    return speciesInfoObj