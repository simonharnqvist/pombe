import pingouin as pg
import pandas as pd


def partial_correlations(df, y, groups):
    """Get the partial correlations between each x and a given y, adjusting for covariates."""
   
    df_numeric = df.apply(pd.to_numeric)
    
    correlations = []
    
    for x in list(df_numeric): # for each variable
    
        # Get which group x belongs to
        x_group = groups.get(x)
        
        # Covariates (z) are the variables that are not x, y, or in the same
        # variable group as x (because that removes most correlations due to
        # expected collinearities)
        z = [col for col in list(df_numeric.drop([x, y], axis=1)) if groups.get(col) != x_group]

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

groups = pd.read_csv("../data/variable_groups.csv")
groups = pd.Series(groups["var_group"].values, index=groups["var"]).to_dict()

correlations = partial_correlations(imputed.drop("Systematic_ID", axis=1), "mean.phylop", groups = groups)
correlations.to_csv("../data/partial_correlations.csv")

