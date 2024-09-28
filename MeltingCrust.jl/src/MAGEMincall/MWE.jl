using MAGEMin_C

data = Initialize_MAGEMin("mp", verbose=false);
Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "MnO"; "H2O"; "O"]
bulk_in_forshaw = [70.48; 25.43/2; 0.77; 3.95; 6.3; 5.54/2; 2.94/2; 0.75/2; 0.07; 100/2; 0.1]
P = 8.0
T = 700.0

sys_in = "mol"

bulk_in_forshaw = [[1], bulk_in_forshaw]
bulk = @view bulk_in_forshaw[2,:]

out     = single_point_minimization(P, T, data, bulk, Xoxides=Xoxides);

grid_model = CreateGrid(nz=100, Lx=2u"km", Lz=28u"km", tfinal=10u"Myr")

@unpack frac_M, bulk, bulk_M, frac_S, bulk_S, ph, ph_frac, ph_type, ph_id, SS_vec, PP_vec, bulk_wt =  out;

save_bulk = copy(bulk)

compo_melt_wt_i = zeros(size(oxides_name,1)-2)

# extract only if melt is more than 7 vol percent
if frac_M >= 0.07
    @unpack Δz, Lx = grid_model
    vol_remove = (100-1/frac_M) / 100
    bulk .= bulk .- vol_remove .* bulk_M .* frac_M

    # store volumed extracted for a 2D surface
    #! check that it works within the function
    magemin_interact.vol_melt_extract_i = frac_M * vol_remove * Lx * Δz

    # renormalize composition without (H2O) and (O)
    index_without_O_H2O =  ("SiO2","Al2O3","CaO","MgO","FeO","K2O","Na2O","TiO2","MnO")

    # find index of O and H2O in out.oxides
    index_O_H2O = indexin(index_without_O_H2O, out.oxides)

    bulk_M_without_H2O = @view bulk_wt[index_O_H2O]

    compo_melt_wt_i .= bulk_M_without_H2O ./ sum(bulk_M_without_H2O)

end

out = Vector{MAGEMin_C.gmin_struct{Float64, Int64}}(undef,2)
data = Initialize_MAGEMin("ig", verbose=false);

function test_MAGEMIN!(out, data)


    P = [10.0, 20.0]
    T = [1100.0, 1200]
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X1 = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    X2 = [49.43; 14.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    X = [X1,X2]
    sys_in = "wt"

    out .= multi_point_minimization(P, T, data, X, Xoxides=Xoxides, sys_in=sys_in);
end

function test_MAGEMIN(data)

    P = [10.0, 20.0]
    T = [1100.0, 1200]
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X1 = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    X2 = [49.43; 14.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    X = [X1,X2]
    sys_in = "wt"

    test = multi_point_minimization(P, T, data, X, Xoxides=Xoxides, sys_in=sys_in);

    return test
end

@btime test_MAGEMIN!(test, data)
# 208.240 ms (1999 allocations: 186.73 KiB)
# 204.927 ms (1531 allocations: 121.23 KiB)

data = Initialize_MAGEMin("ig", verbose=false);
P = [10.0, 20.0]
T = [1100.0, 1200]
Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
X1 = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
X2 = [49.43; 14.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
X = [X1,X2]
sys_in = "wt"

test = Vector{MAGEMin_C.gmin_struct{Float64, Int64}}(undef,2)


@time test .= deepcopy(multi_point_minimization(P, T, data, X, Xoxides=Xoxides, sys_in=sys_in))