import magic
import pandas as pd
from optparse import OptionParser

#Arguments
usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-i", "--csv_input_file", action="store", type="str", dest="csv_input_file",help="csv file of normalized UMIs")
parser.add_option("-o", "--csv_output_file", action="store", type="str", dest="csv_output_file",help="csv output file")
(options, args) = parser.parse_args()

#Run MAGIC
magic_operator = magic.MAGIC()


X = pd.read_csv(options.csv_input_file)
X_magic = magic_operator.fit_transform(X, genes='all_genes')
X_magic.to_csv(options.csv_output_file, index=False)

