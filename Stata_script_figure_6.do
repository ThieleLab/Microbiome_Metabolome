

clear
clear mata
clear matrix
set more off
capture log close

cd P:\AG_Hertel\InSilicoInVivoFiles\Manuscript\Figures

import delimited "P:\AG_Hertel\InSilicoInVivoFiles\Manuscript\Supplement\Table_S1_new.csv"

rename regressioncoefficientsi b_in_silico
rename regressioncoefficientvi b_in_vivo

gen pos=strpos(b_in_silico, "(")-1
gen b_in_silico_per_SD=substr(b_in_silico,1,pos)
replace pos=strpos(b_in_vivo, "(")-1
gen b_in_vivo_per_SD=substr(b_in_vivo,1,pos)
destring(b_in_silico_per_SD b_in_vivo_per_SD), replace
rename name Metabolite
rename fdr FDR

replace Metabolite="Glycerol-3-phosphate" if Metabolite=="Glycerol 3-phosphate"
replace Metabolite="Dodecanoate" if Metabolite=="Dodecanoic acid"
replace Metabolite="Isobutyrate" if Metabolite=="isobutyric acid"
replace Metabolite="Pantothenate" if Metabolite=="Pantothenic acid"
replace Metabolite="Thymidine-5'-phosphate" if Metabolite=="Thymidine 5'-phosphate"

gen met=strupper(substr(Metabolite,1,1))
gen abolite=substr(Metabolite,2,.)
replace Metabolite=met+abolite

local list_metabolites="1,5-diaminopentane 5'-methylthioadenosine Adenine Adenosine Alanine Arginine Asparagine Aspartate Betaine Butyrate Choline Cytidine Dodecanoate Gamma-aminobutyrate Glutamate Glutamine Glutarate Glycerol-3-phosphate Glycine Guanine Guanosine Histamine Histidine Hypoxanthine Inosine Isobutyrate Isoleucine L-Lactate Leucine Lysine Methionine N-acetyl-D-glucosamine Niacinamide Ornithine Pantothenate Phenylalanine Proline Propionate Putrescine Riboflavin Serine Spermidine Succinate Thiamine Threonine Thymidine Thymidine-5'-phosphate Tryptophan Tyramine Tyrosine Urea Uridine Valine"

foreach j in `list_metabolites'{
    local name="`j'"+"_raw_presence"
	twoway (lfitci b_in_vivo_per_SD b_in_silico_per_SD  if Metabolite=="`j'", lcolor(red) clwidth(medthick) clpattern(dash) ciplot(rline) lcolor(red) blwidth(medthin) blpattern(longdash_shortdash) estopts()) (scatter b_in_vivo_per_SD b_in_silico_per_SD if Metabolite=="`j'", mcolor(maroon) msymbol(smdiamond_hollow)) (scatter b_in_vivo_per_SD b_in_silico_per_SD if Metabolite=="`j'" & FDR<0.05, mcolor(navy) msymbol(smdiamond_hollow)) if Metabolite=="`j'", ytitle( "In vivo association statistic" "Δlog c[nmol/g] presence vs absence") ytitle(, size(medsmall)) yline(0, lwidth(medthick) lcolor(black)) ylabel(, nogrid) xtitle("In silico association statistic" "Δj[mmol/d] presence vs absence") xtitle(, size(medium)) xline(0, lwidth(medthick) lcolor(black)) title("`j'", size(vlarge) color(black)) legend(order(1 "95% Confidence interval"  2 "Linear fit (all species)" 3 "Species with FDR>0.05" 4 "Species with FDR<0.05") col(1) region(lcolor(white))) graphregion(fcolor(white) lcolor(white)) saving("`name'", replace) xsize(5.5) ysize(7.5) 
	graph export "`j'_raw_presence.png", replace as(png) 
	}

graph combine isoleucine_raw_presence.gph N-acetyl-D-glucosamine_raw_presence.gph Butyrate_raw_presence.gph Propionate_raw_presence.gph, col(4) xsize(22) ysize(7.5) graphregion(fcolor(white) lcolor(white)) saving(combined_presence_raw, replace) 	
	
clear

import delimited "P:\AG_Hertel\InSilicoInVivoFiles\Manuscript\Supplement\Table_S2_new.csv"	

rename regressioncoefficientsi b_in_silico
rename regressioncoefficientvi b_in_vivo

gen pos=strpos(b_in_silico, "(")-1
gen b_in_silico_per_SD=substr(b_in_silico,1,pos)
replace pos=strpos(b_in_vivo, "(")-1
gen b_in_vivo_per_SD=substr(b_in_vivo,1,pos)
destring(b_in_silico_per_SD b_in_vivo_per_SD), replace
rename name Metabolite
rename fdr FDR

replace Metabolite="Glycerol-3-phosphate" if Metabolite=="Glycerol 3-phosphate"
replace Metabolite="Dodecanoate" if Metabolite=="Dodecanoic acid"
replace Metabolite="Isobutyrate" if Metabolite=="isobutyric acid"
replace Metabolite="Pantothenate" if Metabolite=="Pantothenic acid"
replace Metabolite="Thymidine-5'-phosphate" if Metabolite=="Thymidine 5'-phosphate"

gen met=strupper(substr(Metabolite,1,1))
gen abolite=substr(Metabolite,2,.)
replace Metabolite=met+abolite

cd "P:\AG_Hertel\InSilicoInVivoFiles\Manuscript\Figures\Abundance"	
foreach j in `list_metabolites'{
    local name="`j'"+"_raw_abundance"
	twoway (lfitci b_in_vivo_per_SD b_in_silico_per_SD  if Metabolite=="`j'", lcolor(red) clwidth(medthick) clpattern(dash) ciplot(rline) lcolor(red) blwidth(medthin) blpattern(longdash_shortdash) estopts()) (scatter b_in_vivo_per_SD b_in_silico_per_SD if Metabolite=="`j'", mcolor(maroon) msymbol(smdiamond_hollow)) (scatter b_in_vivo_per_SD b_in_silico_per_SD if Metabolite=="`j'" & FDR<0.05, mcolor(navy) msymbol(smdiamond_hollow)) if Metabolite=="`j'", ytitle( "In vivo association statistic" "Δlog c[nmol/g] per SD abundance") ytitle(, size(medsmall)) yline(0, lwidth(medthick) lcolor(black)) ylabel(, nogrid) xtitle("In silico association statistic" "Δj[mmol/d] per SD abundance") xtitle(, size(medium)) xline(0, lwidth(medthick) lcolor(black)) title("`j'", size(vlarge) color(black)) legend(order(1 "95% Confidence interval"  2 "Linear fit (all species)" 3 "Species with FDR>0.05" 4 "Species with FDR<0.05") col(1) region(lcolor(white))) graphregion(fcolor(white) lcolor(white)) saving("`name'", replace) xsize(5.5) ysize(7.5) 
	graph export "`j'_raw_abundance.png", replace as(png) 
	}
	
graph combine isoleucine_raw_abundance.gph N-acetyl-D-glucosamine_raw_abundance.gph Butyrate_raw_abundance.gph Propionate_raw_abundance.gph, col(4) xsize(22) ysize(7.5) graphregion(fcolor(white) lcolor(white)) saving(abundance_combined_raw, replace) 
	
clear 		