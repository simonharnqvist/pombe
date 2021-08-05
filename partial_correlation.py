import pingouin as pg
import pandas as pd
import argparse
import numpy as np


def partial_correlations(df, y, groups_df):
    """Get the partial correlations between each x and a given y, adjusting for covariates."""
   
    df_numeric = df.apply(pd.to_numeric).replace([np.inf, -np.inf], np.nan).dropna()
    groups_dict = dict(list(zip(groups_df["var"], groups_df["var_group"])))
    correlations = []
    
    for var in groups_df["var"]:
        if var in list(df):
            var_group = groups_dict.get(var)
            z = [n for n in groups_df["var"] if groups_dict.get(n) != var_group]
            z = [n for n in z if n in list(df)]
            if var != y:
                cor = pg.partial_corr(data=df, x=var, y=y, covar=z, method="spearman")
                cor["var"] = var
                cor["group"] = var_group
                correlations.append(cor)
            
    # Combine to single df
    correlations_df = pd.concat(correlations)
            
    # Bonferroni adjust
    correlations_df["p-adj"] = correlations_df["p-val"] * len(correlations)
            
    return correlations_df

parser = argparse.ArgumentParser()
parser.add_argument("--data")
parser.add_argument("--y")
parser.add_argument("--variable_groups")
parser.add_argument("--output")
args = parser.parse_args()

        
def main():
    data = pd.read_csv(args.data).iloc[:,1:].dropna()
    groups = pd.read_csv(args.variable_groups).iloc[:,1:]

    if args.y == "dN":
        data = data.drop(columns = ["mean.phylop"])
    elif args.y == "mean.phylop":
        data = data.drop(columns = ["dN"])
    else:
        raise ValueError("Please specify either dN or mean.phylop as y.") 

    partial_correlations(data.drop(columns = ["Systematic_ID"]), args.y, groups_df = groups).to_csv(args.output, index = None)

if __name__ == "__main__":
    main()


