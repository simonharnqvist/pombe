import pingouin as pg
import pandas as pd

def assign_group(var):
    groups = {"centr":"Network centrality",
                "Function":"Function",
                "Process":"Process",
                "Component":"Component",
                "chr":"Location",
                "CAI":"Codon bias",
                "Residues":"Size",
                "Charge":"Charge",
                "pI":"Charge",
                "Mass..kDa":"Size",
                "sum.protein.cpc":"Expression",
                "sum.mRNA.cpc": "Expression",
                "log.phase.RPKM":"Expression",
                "solid.med.fitness":"Functional importance",
                "essential":"Functional importance",
                "genelength":"Size",
                "start": "Location",
                "end": "Location"}

    # If var len 1 it must be an amino acid
    if len(var) == 1:
        return "AA composition"
    else:
        for key in groups.keys():
            if key in var:
                return groups.get(key)
            else:
                pass




def partial_correlations(df, y):
    """Get the partial correlations between each x and a given y, adjusting for covariates."""
   
    df_numeric = df.apply(pd.to_numeric)
    
    correlations = []
    
    for x in list(df_numeric): # for each variable
    
        # Get which group x belongs to
        x_group = assign_group(x)
        
        # Covariates (z) are the variables that are not x, y, or in the same
        # variable group as x (because that removes most correlations due to
        # expected collinearities)
        z = [col for col in list(df_numeric.drop([x, y], axis=1)) if assign_group(col) != x_group]

        # Get partial correlations
        if x != y:
            cor = pg.partial_corr(data=df_numeric, x=x, y=y,
                                 covar=z, method="spearman")
            cor["var"] = x
            cor["group"] = x_group
            correlations.append(cor)
            
    # Combine to single df
    correlations_df = pd.concat(correlations)
            
    # Bonferroni adjust
    correlations_df["p-adj"] = correlations_df["p-val"] * len(correlations)
            
    return correlations_df

        
# Run and save 
imputed = pd.read_csv("../data/imputed.csv").iloc[:,1:]  
correlations = partial_correlations(imputed.drop("Systematic_ID", axis=1), "mean.phylop")
correlations.to_csv("../data/partial_correlations.csv")


