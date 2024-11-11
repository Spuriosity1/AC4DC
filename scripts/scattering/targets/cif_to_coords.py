#%%
# formatting is patched-together mess, will probably fail often. 


import re




data = """
I1 I 0.21818(4) 0.82462(2) 0.76697(2) 0.04697(12) Uani 1 1 d .
I2 I 0.54054(4) 0.52077(2) 0.87620(2) 0.05084(12) Uani 1 1 d .
I3 I 0.68537(4) 0.86811(2) 0.98562(2) 0.04569(11) Uani 1 1 d .
C8 C 0.4961(5) 0.6505(3) 0.8803(2) 0.0259(9) Uani 1 1 d .
C5 C 0.4560(5) 0.8243(3) 0.8767(2) 0.0242(8) Uani 1 1 d .
C3 C 0.5829(4) 0.7028(3) 0.9214(2) 0.0246(8) Uani 1 1 d .
C7 C 0.3821(4) 0.6830(3) 0.8400(2) 0.0249(8) Uani 1 1 d .
C4 C 0.5620(4) 0.7895(3) 0.9204(2) 0.0238(8) Uani 1 1 d .
C6 C 0.3677(4) 0.7713(3) 0.8369(2) 0.0247(8) Uani 1 1 d .
O1 O 0.6491(4) 0.6365(3) 1.0263(2) 0.0429(9) Uani 1 1 d .
H1 H 0.7093 0.6063 1.0453 0.064 Uiso 1 1 calc R
O3 O 0.3402(4) 0.9557(2) 0.90018(19) 0.0365(8) Uani 1 1 d .
O2 O 0.8231(4) 0.6590(3) 0.9478(2) 0.0521(11) Uani 1 1 d .
O4 O 0.5369(4) 0.9538(2) 0.8328(2) 0.0401(9) Uani 1 1 d .
H4 H 0.5205 1.0044 0.8295 0.060 Uiso 1 1 calc R
N1 N 0.2919(4) 0.6301(3) 0.8012(2) 0.0365(9) Uani 1 1 d .
H1A H 0.3057 0.5764 0.8017 0.044 Uiso 1 1 calc R
H1B H 0.2227 0.6511 0.7769 0.044 Uiso 1 1 calc R
C1 C 0.6990(5) 0.6635(3) 0.9665(3) 0.0310(10) Uani 1 1 d .
C2 C 0.4382(5) 0.9184(3) 0.8717(2) 0.0265(9) Uani 1 1 d .
O1W O 0.5363(4) 0.1151(2) 0.8245(2) 0.0425(9) Uani 1 1 d D
H1WA H 0.5977 0.1194 0.7912 0.051 Uiso 1 1 d RD
H1WB H 0.5755 0.1315 0.8627 0.051 Uiso 1 1 d RD
"""

a,b,c = 9.3044, 15.8623, 18.9924


output = []
atom_id = 1
ligand = "I3C"

for line in data.strip().splitlines():
    parts = re.split(r'\s+', line)
    atom_name = parts[0]
    atom_type = parts[1]
    x = float(re.sub(r'\(.*\)', '', parts[2]))*a
    y = float(re.sub(r'\(.*\)', '', parts[3]))*b
    z = float(re.sub(r'\(.*\)', '', parts[4]))*c

    output_line = f"ATOM{atom_id:{' '}>7}  {atom_name:<4}{ligand} A 500     {x:>7.3f} {y:>7.3f} {z:>7.3f}  1.00  0.00"
    output.append(output_line)
    atom_id += 1


result = "\n".join(output)
print(result)
# %%
