#%%
import pandas as pd
from pyliftover import LiftOver
df = pd.read_excel("ep_enhancer_atlas_2.0_mcf-7.xlsx")
lo = LiftOver('hg19ToHg38.over.chain')
# %%
def liftover_coordinates(chrom,coord):
    lift_coord = lo.convert_coordinate(chrom, coord)
    if lift_coord is not None and len(lift_coord) == 1:
        return lift_coord[0][0], lift_coord[0][1]
    else:
        return None, None
# %%
# Apply the lift_coordinates function to the DataFrame
df[["lifted_e_chr_start", "lifted_e_start"]] =  df.apply(lambda row: pd.Series(liftover_coordinates(row['e_chr'], row['e_start'])), axis=1)
df[["lifted_e_chr_end", "lifted_e_end"]] =  df.apply(lambda row: pd.Series(liftover_coordinates(row['e_chr'], row['e_end'])), axis=1)
df[["lifted_gene_chr", "lifted_gene_tss"]] =  df.apply(lambda row: pd.Series(liftover_coordinates(row['gene_chr'], row['gene_tss'])), axis=1)

# %%
df.to_excel("lifted_ep_mcf7_enhanceratlas_pyliftover.xlsx", header = True, index = False)
# %%

# %%
