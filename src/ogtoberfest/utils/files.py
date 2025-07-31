import os
import pathlib
from ogtoberfest.utils import util



class FileHandler:
    def __init__(
            self, 
            database_path: pathlib.Path,
            working_dir: pathlib.Path,
            q_dna: bool = False
            # method: str = "OrthoFinder3",
            # msa_method: str = "mafft",
            # tree_method: str = "fasttree",
    ):
        self.wd_base = working_dir
        self.fasta_dir = database_path
        self.q_dna = q_dna

        # self.wd_current = working_dir / method

        # if not os.path.exists(self.wd_current):
        #     os.makedirs(self.wd_current, exist_ok=True)

        # self.msa_dir = self.wd_current / msa_method 
        # if not os.path.exists(self.msa_dir):
        #     os.makedirs(self.msa_dir, exist_ok=True)

        # self.tree_dir = self.wd_current / tree_method
        # if not os.path.exists(self.tree_dir):
        #     os.makedirs(self.tree_dir, exist_ok=True)  


    def get_species_ids_fn(self): 
        return self.wd_base / "SpeciesIDs.txt"
        
    def get_sequence_ids_fn(self):
        # It is always in the first of the 'extension' directories (as this is the relevant one)
        return self.wd_base / "SequenceIDs.txt"

    # def get_msa_dir(self):
    #     return self.msa_dir 
    
    # def get_tree_dir(self):
    #     return self.tree_dir

    # def get_fasta_dir(self):
    #     return self.fasta_dir

    def get_species_fasta_fn(self, iSpecies, qForCreation=False):
        """
        qForCreation: A path is required at which the file should be created (don't search for it)
        """

        fn = self.wd_base / f"Species{iSpecies}.fa" 
        if qForCreation:
            return fn
        if os.path.exists(fn): 
            return fn
        raise Exception(fn + " not found")

    
    def process_fasta_files(self):

        speciesInfoObj = util.SpeciesInfo()

        fasta_extensions = {"fa", "faa", "fasta", "fas", "pep", "fna"}
        qOk = True
        if not os.path.exists(self.fasta_dir):
            print(f"\nDirectory does not exist: {self.fasta_dir}")
            util.Fail()

        files_in_directory = sorted([
            f for f in os.listdir(self.fasta_dir) if os.path.isfile(os.path.join(self.fasta_dir, f))
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
            print("OrthoFinder expects FASTA files to have one of the following extensions: %s" % (", ".join(fasta_extensions)))
        
        if len(fasta_files) < 2:
            print("ERROR: At least two species are required")
            if len(fasta_files) == 0:
                print(f"\nNo fasta files found in supplied directory: {self.fasta_dir}")
            util.Fail()

        iSeq = 0
        iSpecies = 0

        with open(self.get_sequence_ids_fn(), 'a') as sequence_ids_file, \
        open(self.get_species_ids_fn(), 'a') as species_ids_file:
            for fastaFilename in fasta_files:
                speciesInfoObj.speciesToUse.append(iSpecies)

                outputFasta = open(self.get_species_fasta_fn(iSpecies, qForCreation=True), 'w')
                fastaFilename = fastaFilename.rstrip()
                species_ids_file.write("%d: %s\n" % (iSpecies, fastaFilename))
                baseFilename, extension = os.path.splitext(fastaFilename)
                speciesInfoObj.speciesToUseName.append(baseFilename)

                species_name = (
                    baseFilename.split(".")[0].replace(":", "_")
                    .replace(",", "_")
                    .replace("(", "_")
                    .replace(")", "_")
                    .replace(" ", "_")
                )
                speciesInfoObj.species2id_dict[species_name] = iSpecies
                speciesInfoObj.id2species_dict[iSpecies] = species_name
                
                mLinesToCheck = 100
                qHasAA = False
                with open(self.fasta_dir / fastaFilename, 'r') as fastaFile:
                    for iLine, line in enumerate(fastaFile):
                        if line.isspace(): continue
                        if len(line) > 0 and line[0] == ">":
                            newID = "%d_%d" % (iSpecies, iSeq)
                            acc = line[1:].rstrip()
                            if len(acc) == 0:
                                print("ERROR: %s contains a blank accession line on line %d" % (self.fasta_dir / fastaFilename, iLine+1))
                                util.Fail()
                            sequence_ids_file.write("%s: %s\n" % (newID, acc))

                            sequence_id = newID.replace("#", "").strip()

                            accession = (
                                acc.replace(":", "_")
                                .replace(",", "_")
                                .replace("(", "_")
                                .replace(")", "_")
                            )

                            # if accession.split(".")[0].lower() != baseFilename.split(".")[0].lower():
                            #     accession = baseFilename.split(".")[0] + "." + accession
                            accession = accession.split()[0]
                            if species_name.lower() not in accession.lower():
                               accession = species_name + "." + accession

                            speciesInfoObj.sequence2id_dict[accession] = sequence_id
                            speciesInfoObj.id2sequence_dict[sequence_id] = accession

                            outputFasta.write(">%s\n" % newID)    
                            iSeq += 1
                        else:
                            line = line.upper()    # allow lowercase letters in sequences
                            if not qHasAA and (iLine < mLinesToCheck):
    #                            qHasAA = qHasAA or any([c in line for c in ['D','E','F','H','I','K','L','M','N','P','Q','R','S','V','W','Y']])
                                qHasAA = qHasAA or any([c in line for c in ['E','F','I','L','P','Q']]) # AAs minus nucleotide ambiguity codes
                            outputFasta.write(line)
                    outputFasta.write("\n")
                if (not qHasAA) and (not self.q_dna):
                    qOk = False
                    print("ERROR: %s appears to contain nucleotide sequences instead of amino acid sequences. Use '-d' option" % fastaFilename)
                iSpecies += 1
                iSeq = 0
                outputFasta.close()
            if not qOk:
                util.Fail()

        if len(fasta_files) > 0: 
            outputFasta.close()

        speciesInfoObj.nSpAll = max(speciesInfoObj.speciesToUse) + 1      # will be one of the new species
    
        return speciesInfoObj



