
import tempfile
import subprocess

class CDHitGroup:

    #./cd-hit -i input.txt -o output.txt -c 0.9 -n 5 -M 1000 -d 0 -T 1 -bak 1
    def get_group(self, proteins, similarity=0.9, word_size=5, memory_mb=1024, threads=1):
        """Given a list of pyprot.Protein objects, runs CDHit
        to cluster sequences based on similarity.
        CDHit must be installed for this method to work."""
        # what's up with the sequence of a protein with many chains?
        sequences = [prot.get_sequence() for prot in proteins]
        #write in a tempdir file hopefully.
        # use input.txt and etc in tempdir
        # check where cdhit is
        # assert that CDHit is installed..
        out = subprocess.run(["./cd-hit",
            "-i input.txt", "-o output.txt",
            "-c {}".format(similarity), "-n {}".format(word_size),
            "-M {}".format(memory_mb), "-T {}".format(threads),
            "-bak 1", "-d 0"], capture_output=True)
