#!/usr/bin/bash
base_path=`git rev-parse --show-toplevel`
echo ${base_path}

# Download expression data
wget ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/annotation/TAR_gene_expression/c_elegans.PRJNA13758.current.TAR_gene_expression.tar.gz -O ${base_path}/data/expression.gz
extract ${base_path}/data/expression.gz

