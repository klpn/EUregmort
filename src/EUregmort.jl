module EUregmort

using CSV, DataFrames, PyCall, PyPlot, Statistics

const shpreader = PyNULL()
const ccrs = PyNULL()

function __init__()
	copy!(shpreader, pyimport("cartopy.io.shapereader"))
	copy!(ccrs, pyimport("cartopy.crs"))
end

export nuts2ids, caprop_regcmp, caprop_regsexplot, fourp, fivep, caprop_mapplot, meanrate
mainpath = normpath(@__DIR__, "..")
datapath = normpath(mainpath, "data")
shpdatapath = normpath(datapath, "NUTS_2013_03M_SH", "data") 
ycdr = CSV.read(normpath(datapath, "hlth_cd_ycdr2.csv"); missingstring = ":", rows_for_type_detect = 200)
nuts = CSV.read(normpath(datapath, "NUTS_AT_2013.csv"); rows_for_type_detect = 200)
nids = convert(Array, collect(skipmissing(nuts[:NUTS_ID])))
ageseq = ["Y_LT1"; "Y1-4"; "Y5-9"; "Y10-14"; "Y15-19"; 
	"Y20-24"; "Y25-29"; "Y30-34"; "Y35-39"; "Y40-44"; "Y45-49"; 
	"Y50-54"; "Y55-59"; "Y60-64"; "Y65-69"; "Y70-74"; "Y75-79"; 
	"Y80-84"; "Y85-89"; "Y90-94"; "Y_GE95"]
agesplitter(age) = split(age,  ['Y'; 'M'; '_'; '-'], keepempty = false)

function agealias(age)
	aspl = agesplitter(age)
	if aspl[1] == "TOTAL"
		alias = "alla åldrar"
	elseif aspl[1] == "LT1"
		alias = "0"
	elseif startswith(aspl[1], "LT")
		alias = "0\u2013$(parse(aspl[end][3:end])-1)"
	elseif startswith(aspl[1], "GE")
		alias = "$(aspl[end][3:end])\u2013"
	else
		alias = "$(aspl[1])\u2013$(aspl[end])"
	end

	if (startswith(age, "YM") || startswith(age, "M"))
		return "$alias genomsnitt över åldrar"
	else
		return alias
	end
end

sexaliases = Dict("F" => "kvinnor", "M" => "män", "T" => "alla")
causealiases = Dict("pop" => L"$10^5$ personår", "A-R_V-Y" => "alla orsaker")

function causealias(ca)
	if ca in keys(causealiases)
		return causealiases[ca]
	else
		return ca
	end
end

function nuts2ids(ctry) 
	n2ids = map((x)->Symbol(x), filter((x)->(length(x)==4 && startswith(x, ctry)), nids))
	filter!((x) -> x in names(ycdr), n2ids)
end

dfarrmatch(col, arr) = map((x) -> in(x, arr), Vector(col))


function meanrate(sage, eage)
	sexframes = Dict()
	sind = indexin([sage], ageseq)[1]
	eind = indexin([eage], ageseq)[1]
	agelist = ageseq[sind:eind]
	if sage == "Y_LT1"
		if eage == "Y_GE95"
			agestr = "MTOTAL"
		else
			agestr = "YM_LT$(parse(agesplitter(eage)[end])+1)"
		end
	else
		if eage == "Y_GE95"
			agestr = "YM_GE$(agesplitter(sage)[1])"
		else
			agestr = "YM$(agesplitter(sage)[1])-$(agesplitter(eage)[end])"
		end
	end
	for sex in keys(sexaliases)
		ycdr_sub = ycdr[((ycdr[:sex].==sex) & (dfarrmatch(ycdr[:age], agelist))),
			 [4;6:size(ycdr)[2]]]
		ycdr_sub_agg = aggregate(ycdr_sub, :icd10, mean)
		aggregs = rename((x)->Symbol(replace(string(x), "_mean", "")),
			ycdr_sub_agg[:, 2:end])
		aggl = size(ycdr_sub_agg)[1]
		metaframe = DataFrame(unit = fill("RT", aggl), sex = fill(sex, aggl),
			age = fill(agestr, aggl),
			icd10 = ycdr_sub_agg[:icd10], time_geo = fill(2013, aggl)) 
		sexframes[sex] = hcat(metaframe, aggregs) 
	end
	vcat(values(sexframes)...)
end

function ycdr_regcmp(nutsregs, sex, age, icd10, time_geo, inframe)
	inframe_sub = DataFrame(unit = inframe[:unit], sex = inframe[:sex],
		age = inframe[:age], icd10 = inframe[:icd10],
		time_geo = inframe[:time_geo])
	for reg in nutsregs
		inframe_sub[reg] = inframe[reg]
	end
	inframe_long = stack(inframe_sub, nutsregs)
	inframe_long[((inframe_long[:sex].==sex) .& (inframe_long[:age].==age)
		.& (inframe_long[:icd10].==icd10) .& (inframe_long[:time_geo].==time_geo)), :]
end

function caprop_regcmp(nutsregs, sex, age, ca1, ca2, time_geo, inframe)
	ca1frame = ycdr_regcmp(nutsregs, sex, age, ca1, time_geo, inframe)
	if ca2 == "pop"
		capropframe = ca1frame
	else
		ca2frame = ycdr_regcmp(nutsregs, sex, age, ca2, time_geo, inframe)
		capropframe = copy(ca1frame)
		capropframe[:value] = ca1frame[:value]./ca2frame[:value]
	end
	capropframe
end

function caprop_regsexplot(nutsregs, age, ca1; ca2 = "A-R_V-Y", time_geo = 2013,
	inframe = ycdr)
	fframe = caprop_regcmp(nutsregs, "F", age, ca1, ca2, time_geo, inframe)
	mframe = caprop_regcmp(nutsregs, "M", age, ca1, ca2, time_geo, inframe)
	scatter(fframe[:value], mframe[:value], alpha=0.3)
	for i in 1:size(nutsregs)[1]
		text(fframe[:value][i], mframe[:value][i], fframe[:variable][i],
			ha = "center", va = "center")
	end
	grid(1)
	xlabel(sexaliases["F"])
	ylabel(sexaliases["M"])
	title(*("Döda $(time_geo) $(causealias(ca1))/$(causealias(ca2)) ",
		"$(agealias(age))"))
end

perc_round(value) = replace("$(round(value; digits=4))", "." => ",")
fourp(prop) = 
	[
	Dict("col" => "lightyellow", "value" => quantile(prop, 1/4));
	Dict("col" => "yellow", "value" => quantile(prop, 2/4));
	Dict("col" => "tomato", "value" => quantile(prop, 3/4));
	Dict("col" => "red", "value" => quantile(prop, 1))
	]
fivep(prop) = 
	[
	Dict("col" => "lightyellow", "value" => quantile(prop, 1/5));
	Dict("col" => "yellow", "value" => quantile(prop, 2/5));
	Dict("col" => "orange", "value" => quantile(prop, 3/5));
	Dict("col" => "tomato", "value" => quantile(prop, 4/5));
	Dict("col" => "red", "value" => quantile(prop, 1))
	]

function caprop_mapplot(nutsregs, sex, age, ca1; ca2 = "A-R_V-Y", time_geo = 2013,
		percfunc = fivep,
		shapefname = normpath(shpdatapath, "NUTS_RG_03M_2013_3034.shp"),
		inframe = ycdr)
	pframe = caprop_regcmp(nutsregs, sex, age, ca1, ca2, time_geo, inframe)
	region_shp = shpreader[:Reader](shapefname)
	regstrings = map((x)->string(x), nutsregs)
	proj = ccrs[:LambertConformal](central_longitude = 10, central_latitude = 52,
		standard_parallels = (35,65), false_easting = 4000000,
		false_northing = 2800000, globe = ccrs[:Globe](ellipse = "GRS80"))
	ax = plt[:axes](projection = proj)
	prop = pframe[:value]
	propdict = Dict(zip(regstrings, prop))
	quantiles = percfunc(prop) 
	boundlist = []
	fcolor = "blue"
	for region_rec in region_shp[:records]()
		nid = region_rec[:attributes]["NUTS_ID"]  
		if nid in regstrings
			boundlist = vcat(boundlist, region_rec[:bounds])
			xmean = mean([boundlist[end][1];
				boundlist[end][3]])
			ymean = mean([boundlist[end][2];
				boundlist[end][4]])
			for quantile in quantiles
				if propdict[nid] <= quantile["value"]
					fcolor = quantile["col"]
					break
				end
			end
			ax[:add_geometries](region_rec[:geometry], proj,
				edgecolor = "black", facecolor = fcolor)
			ax[:annotate](nid, (xmean, ymean), ha = "center")
		end
	end
	xminimum = minimum([bound[1] for bound in boundlist])
	xmaximum = maximum([bound[3] for bound in boundlist])
	yminimum = minimum([bound[2] for bound in boundlist])
	ymaximum = maximum([bound[4] for bound in boundlist])
	ax[:set_xlim](xminimum, xmaximum)
	ax[:set_ylim](yminimum, ymaximum)
	percpatches = []
	perclabels = []
	for (i, quantile) in enumerate(quantiles)
		percpatch = matplotlib[:patches][:Rectangle]((0, 0), 1, 1,
			facecolor = quantile["col"])
		percpatches = vcat(percpatches, percpatch)
		if i == 1
			perclabel = *("\u2265", perc_round(minimum(prop)),
				"\n\u2264", perc_round(quantile["value"]))
		else
			perclabel = *("\u2264", perc_round(quantile["value"]))
		end
		perclabels = vcat(perclabels, perclabel)
	end
	legend(percpatches, perclabels, loc = "upper left", 
		framealpha = 0.75, bbox_to_anchor=(1,1))
	title(*("Döda $(time_geo) $(causealias(ca1))/$(causealias(ca2)) ",
		"$(sexaliases[sex]) $(agealias(age))"))
	subplots_adjust(right = 0.8)
end

end # module
