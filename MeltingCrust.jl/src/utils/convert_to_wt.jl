# using Parameters

# convert mole of element in wt of oxide
# function convert_compo_to_wt_dry_pct!(compo_i_wt, compo_i_n)

#     # define gram formula weight of each oxide
#     SiO2_GFW = 60.0843
#     Al2O3_GFW = 101.9613
#     CaO_GFW = 56.0774
#     MgO_GFW = 40.3044
#     FeO_GFW = 71.8444
#     K2O_GFW = 94.1960
#     Na2O_GFW = 61.9789
#     TiO2_GFW = 79.8664
#     MnO_GFW = 70.9374
#     # H2O_GFW = 18.0153
#     # GFW_with_ratio = (SiO2_GFW, Al2O3_GFW/2, CaO_GFW, MgO_GFW, FeO_GFW, K2O_GFW/2, Na2O_GFW/2, TiO2_GFW, MnO_GFW, H2O_GFW/2)
#     # without H2O because TAS is dry
#     GFW_with_ratio = (SiO2_GFW, Al2O3_GFW, CaO_GFW, MgO_GFW, FeO_GFW, K2O_GFW, Na2O_GFW, TiO2_GFW, MnO_GFW)

#     compo_i_n_wo_O = @view compo_i_n[1:end-2]

#     compo_i_wt .= compo_i_n_wo_O .* GFW_with_ratio ./ sum(compo_i_g)
# end


function convert_compo_to_wt_dry!(compo_i_wt_dry, compo_i_wt)
    # renormalize compo_i_wt by removing 2 last indexes
    compo_i_wt_wo_O = @view compo_i_wt[1:end-2]
    compo_i_wt_dry .= compo_i_wt_wo_O ./ sum(compo_i_wt_wo_O)
end



# # convert mole of element in wt of oxide
function convert_compo_to_wt_pct!(compo_i_wt, compo_i_n, compo_i_g)

    # define gram formula weight of each oxide
    SiO2_GFW = 60.0843
    Al2O3_GFW = 101.9613
    CaO_GFW = 56.0774
    MgO_GFW = 40.3044
    FeO_GFW = 71.8444
    K2O_GFW = 94.1960
    Na2O_GFW = 61.9789
    TiO2_GFW = 79.8664
    MnO_GFW = 70.9374
    H2O_GFW = 18.0153
    O_GFW = 16.0
    GFW_with_ratio = (SiO2_GFW, Al2O3_GFW, CaO_GFW, MgO_GFW, FeO_GFW, K2O_GFW, Na2O_GFW, TiO2_GFW, MnO_GFW, H2O_GFW, O_GFW)

    compo_i_g .= compo_i_n .* GFW_with_ratio
    compo_i_wt .= compo_i_g ./ sum(compo_i_g)
end


function convert_compo_solid_to_wt_pct!(compo_rock_wt, compo_solid_matrix, compo_i_g)

    # convert compo_rock_wt from mol% to wt%
    for i in axes(compo_solid_matrix,1)
        compo_rock_n_i = @view compo_solid_matrix[i,:]
        compo_rock_wt_i = @view compo_rock_wt[i,:]
        convert_compo_to_wt_pct!(compo_rock_wt_i, compo_rock_n_i, compo_i_g)
    end
end