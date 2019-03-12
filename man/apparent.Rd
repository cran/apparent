\name{apparent}
\alias{apparent}
\alias{Parentage analysis}
\title{
	Accurate parentage analysis in the absence of guiding information
}
\description{ 
Performs parentage analysis based on a test of genetic identity between expected progeny (EPij), built using Single Nucleotide Polymorphism (SNP) homozygous loci from all pairs of possible parents (i and j), and all potential offspring (POk). Using the Gower Dissimilarity metric (GD), genetic identity between EPij and POk is taken as evidence that individuals i and j are the true parents of offspring k. Evaluation of triad (two parents + offspring) significance is based on the distribution of all GD (EPij|POk) values. Specifically, a Dixon test is used to identify a gap-based threshold that separates true triads and from spurious associations. For any offspring not successfully assigned to a pair of parents, perhaps due to the absence of one parent from the test population, a non-mandatory Dyad analysis can be employed to identify a likely single parent for a given offspring. In this analysis, a two-stage test is applied to discriminate an offspring's true parent from its other close relatives (e.g. siblings) that may also be present in the population. In the first stage, 'apparent' calculates the mean GD (GDM) between a POk and all expected progeny arising from the j possible triads involving potential parent i. In the second stage, it calculates a coefficient of variation (GDCV) among the pairwise GD's between POk and each expected progeny arising from the j triads involving potential parent i. An individual that is simultaneously a low outlier in the first test and a high outlier in the second is identified as a likely parent of POk. In an effort to facilitate interpretation, results of both the triad and optional dyad analyses are presented in tabular and graphical form.
}
\usage{
apparent(InputFile, MaxIdent=0.10, alpha=0.01, nloci=300, 
self=TRUE, plot=TRUE, files=TRUE, Dyad=FALSE, dir)
}
\arguments{
  \item{InputFile}{
A tab-delimited text file with n rows (n = number of individuals in the population) and m+2 columns (m = number of SNP loci). Column one contains the ID's of the individuals in the population.  Column two contains a classification key which assigns each individual to one of five possible classes from analysis (see below). The third and all subsequent columns contain genotype calls, with one SNP per column and the alleles separated by "/".  A missing genotype is represented by "-/-".  An example input file with 5 individuals and 5 SNP loci is shown below:
\cr
\cr     Geno1   All   A/A   A/T   C/A   C/C   T/C
\cr     Geno2   Pa    T/A   A/A   -/-   C/C   T/T
\cr     Geno3   Mo    A/A   A/A   C/A   C/G   T/T
\cr     Geno4   Fa    T/T   A/T   A/A   G/G   -/-
\cr     Geno5   Off   A/A   T/T   C/C   C/C   T/C
\cr
\cr 	The keys values allowed in the second column are: 
\cr		All: "All" - The individual will be tested as a potential Mother, Father, and Offspring.
\cr		Pa: "Parent" - The individual will be tested as potential Mother and Father.
\cr		Mo: "Mother" - The individual will be tested only as a potential Mother.
\cr		Fa: "Father" - The individual will be tested only as a potential Father.
\cr		Off: "Offspring" - The individual will be tested only as a potential Offspring.
  }
  \item{MaxIdent}{
Sets the maximum triad GDij|k to be considered for outlier significance testing.  This parameter directly impacts computation time.  By default, MaxIdent is set to 0.1. 
  }
  \item{alpha}{
The alpha level for all significance testing (triad and dyad analyses).
  }
  \item{nloci}{
The minimum acceptable number loci to be used when computing the pairwise GDij|k.  The default value of 300 is suggested, based on previous investigations.  All triads for which the number of usable SNPs falls below nloci will be excluded from the analysis.
  }
  \item{self}{
Logical value for instructing 'apparent' whether or not to consider self-crossing (parent i = parent j). The default value is TRUE.
  }
  \item{plot}{
Logical value for create the plots of both Triad and Dyad analysis. By default, 'apparent' creates only the Triad analysis plot. The default value is TRUE, meaning 'apparent' creates only the Triad analysis plot. The Dyad analysis plot can't be printed, only saved (it is a pdf file with multiple plots for each significant offspring). To save Dyad analysis plot, make sure plot, files and Dyad flags are TRUE.
  }
  \item{files}{
Logical value for saving the tables and graphics output files from both Triad and Dyad analysis. The default value is FALSE, meaning the results will be stored at local variables, instead saved. The Triad analysis files that can be saved are: apparent-Triad-All.csv, apparent-Triad-Sig.csv, apparent-Triad-Summary.csv and apparent-Triad-Plot.png. The Dyad analysis result table (apparent-Dyad-Sig.csv) will be saved only if both files and Dyad flags are TRUE. To save Dyad analysis plot, make sure plot, files and Dyad flags are TRUE.  
  }
  \item{Dyad}{
Logical value for instructing 'apparent' to perform a Dyad analysis, following the Triad analysis. The default value is FALSE.  
  }
  \item{dir}{
A path to work directory where the files should be saved. MacOS example for Desktop savings: "/Users/username/Desktop.  
  }
}
\details{
See Melo ATO and Hale I (2019) 'apparent': a simple and flexible R package for accurate SNP-based parentage analysis in the absence of guiding information. BMC Bioinformatics. doi: doi.org/10.1186/s12859-019-2662-3. Also, find useful information at https://github.com/halelab/apparent.
}
\note{
As the size of the exploratory triad space inflates (i.e. as the population size grows), the likelihood of low-probability instances of low GD values, occurring purely by chance, increases.  Such outlier values can "bleed" into the gap between true and spurious triads, in effect making the gap disappear and leading to a null result, even if true triads exist in the population (Type II error). Such a dilution of a true gap is generally due to the presence of one or more problematic genotypes in the population, problematic due to the fact that the available number of loci for their analysis is much for some reason lower than for other individuals, leading to artificially low GDij|POk values. To address such problems the following two "best practices" are recommended to users: 1. Always include at least one known (true) triad in the analysis, as a reference/control. 2. After running the triad analysis, inspect the summary file "apparent-Triad-Summary.cvs" to see if any individuals exhibit unusually low mean GDij|k values or mean number of loci usable for analysis. If such individuals exist, remove them and re-run the analysis.
}
\value{
"Triad_All" and ""apparent-Triad-All.csv" contains the full table of results from the triad analysis, including the cross type (self or out-cross), the number of usable SNPs, and the GDij|k, of all tested triads.

"Triad_Log1", "TriadLog2" and ""apparent-Triad-Summary.csv" reports useful summary statistics from the triad analysis. For the population as a whole: Overall mean GDij|k and its standard deviation, as well as overall mean number of usable SNPs and its standard deviation. For each individual genotype in the analysis: Mean GDij|k, GDij|k range, and mean usable number of SNPs for all comparisons involving that genotype. Genotypes exhibiting a mean GDij|k or number of usable SNPs less than 2 SD's below the population means are considered outliers.

"Triad_Sig" and "apparent-Triad-Sig.csv" reports only those triads found to be significant.

"apparent-Triad-Plot.png" is a plot of the distribution of GDij|k values, annotated with the gap-based threshold that separates true triads from spurious associations.  The plot is useful in interpretating the results of the triad analysis , the full details of which are found in "apparent-Triad-All.txt" and "apparent-Triad-Sig.txt".

"Dyad_Sig" and "apparent-Dyad-Sig.csv" reports only those parent-offspring dyads found to be significant.

"apparent-Dyad-Plot.pdf" is a two-paneled figure showing the distributions of GDM and GDCV values upon which the dyad analysis is based.  These plots are useful in interpretating the results of the dyad analysis, the full details of which are found in "apparent-Dyad-Sig.txt".
}
\references{ 
Dixon W.J. (1950) Analysis of extreme values. Ann. Math. Stat. 21(4):488-506.
Dixon W.J. (1951) Ratios involving extreme values. Ann. Math. Stat. 22(1):68-78.
Gower JC. (1971) A General Coefficient of Similarity and Some of Its Properties. Biometrics. 27(4):857-71.
}
\author{
 Arthur T. O. Melo and Iago Hale
Department of Agriculture Nutrition and Food Systems, University of New Hampshire, Durham, NH, USA.
}
\seealso{ \code{\link{outlier}}}
\examples{
# Load the input file
InputFile <- apparent_TestData
# Run the apparent
apparent(InputFile,dir=tempdir())
}
\keyword{Parentage analysis}
\keyword{Pedigree inference}
\keyword{R package}
\keyword{Single Nucleotide Polymorphism (SNP)}