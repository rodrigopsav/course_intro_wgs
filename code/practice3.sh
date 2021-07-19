#-------------------------------------#
#          Practice 3: IVDP           #
#-------------------------------------#


######## Running IVDP ########

### Run on a local server ####
# Example 2
./ivdp.sh -p params/pe_ex2.txt

# Example 3b
./ivdp.sh -p params/pe_ex3b.txt

# Example 10
./ivdp.sh -p params/se_ex10.txt

# Example 12
./ivdp.sh -p params/se_ex12.txt


### Run on HPCC with slurm ####
# Example 2
./ivdp.sh -p params/pe_ex2.txt -c params/configSlurm.txt

# Example 3b
./ivdp.sh -p params/pe_ex3b.txt -c params/configSlurm.txt

# Example 10
./ivdp.sh -p params/se_ex10.txt -c params/configSlurm.txt

# Example 12
./ivdp.sh -p params/se_ex12.txt -c params/configSlurm.txt


# Note: To open html and pdf files with command line, use:
# evince file.pdf
# firefox file.html
