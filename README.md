# DriverMP

In this study, we present a new method termed DriverMP (Multiomics-based Pair driver genes) for effectively prioritizing altered genes on a cancer type level by considering mutation pair genes. It is designed by first applying somatic mutation data, protein-protein interaction network, and differential gene expression data to prioritizing mutation pairs, based on which individual mutated genes are then prioritized.

	
## Requirements for installation

**A.** `Linux/Unix` environment


**B.** `C++` is one of the program languages, which is available for most operating systems. If you are in UNIX/Linux, it is probably already installed;


**C.** We tested the code on `Ubuntu 18.04.6 LTS`, using `g++ version 7.5.0`, and we recommend that you use DriverMP on that version or above. If you want to check your current version of g++, you can use the following command in your terminal:


`g++ --version`


**D.** Unzip "Driver.zip" you have downloaded;


**E.** You can use our compiled `DriverMP` directly as shown below, 


or you can compile it yourself using the code we provide as follows：


`g++ -Ofast -fopenmp *.cpp -o DriverMP`


Then, you can use DriverMP as normal.

## Running DriverMP
		
**A.** You can run DriverMP as follow.

`./DriverMP [Mutation data] [Tumor expression data] [Normal expression data] [PPI network data]`


*[Mutation data]: Path to format-compliant non-silent somatic mutation data;*

*[Tumor expression data]: Path to RNA-Seq data (FPKM normalized) for specific cancer samples;*

*[Normal expression data]: Path to RNA-Seq data (FPKM normalized) of normal samples corresponding to the specific cancer;*

*[PPI network data]: Path to Protein-Protein Interaction Networks;*


**B.** You can try the default example we have provided as follows.

`./DriverMP`

*This example will run directly on the breast cancer data and HumanNet interaction network data we provided.*

				
**C.** You can provide only cancer-specific data and use our default PPI network.
				
`./DriverMP [Mutation data] [Tumor expression data] [Normal expression data] H`

or

`./DriverMP [Mutation data] [Tumor expression data] [Normal expression data] S`

*H and S represent the HumanNet and STRINGv10 networks that we provide for you by default, respectively.*
		

**D.** Output

The output results are saved to the `Output//output.txt` file by default.
	
## Data Format Requirements


**A.** Mutation Data

(i) Note that only genes with `Official Symbol` can be recognised;

(ii) The first row is the sample ID, and the first element of each row starting from the second row is the gene of form Official Symbol.

The format requirement is as follows:

**[sample 1]** tab **[sample 2]** tab … tab **[sample M]**

**[  Gene  A ]** tab [ 0 or 1 ] tab [ 0 or 1 ] tab … tab [ 0 or 1 ]

**[  Gene  B ]** tab [ 0 or 1 ] tab [ 0 or 1 ] tab … tab [ 0 or 1 ]

...
	    
**[  Gene  C ]** tab [ 0 or 1 ] tab [ 0 or 1 ] tab … tab [ 0 or 1 ]


**B.** Gene Expression RNA-Seq Data

(i) Note that only genes with `Official Symbol` can be recognised;

(ii) The first row is the sample ID, and the first element of each row starting from the second row is the gene of form Official Symbol;

(iii) We use RNA-Seq data downloaded in FPKM format in TCGA, and we recommend that you use this format as well;

(iv) You need to provide both `Tumor and normal` RNA-Seq data, please make sure they both conform to the format.

The format requirement is as follows:

**[sample 1]** tab **[sample 2]** tab … tab **[sample M]**

**[  Gene  A ]** tab [ value ] tab [ value ] tab … tab [ value ]

**[  Gene  B ]** tab [ value ] tab [ value ] tab … tab [ value ]

...
	    
**[  Gene  C ]** tab [ value ] tab [ value ] tab … tab [ value ]


**C.** Protein-Protein Interaction Network

(i) DriverMP provides two default PPI networks HumanNet and STRINGv10, if you need to use other PPI networks please provide formatted data.

(ii) At present, the node naming method of PPI network is generally in NCBI format. For your convenience, we have built-in conversion of genes' names, and you can directly provide the downloaded PPI network to DriverMP.

(iii) For the best performance of DriverMP, you should preferably provide PPI networks with weights that are maximum-minimum normalized as follows: 

$$
x^{'} = \frac{x - \min x}{\max x - \min x}
$$
		
The data format requirements are as follow: 
			
**[Gene NCBI A]** tab **[Gene NCBI ID B]** tab [the weight between Gene A and B]

...

**[Gene NCBI C]** tab **[Gene NCBI D]** tab [the weight between Gene C and D]
