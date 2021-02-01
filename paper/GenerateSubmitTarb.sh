#! /bin/bash
# script to generate submission tarball
# it creates a dir with tex and bib source and fig
# figs are renamed and numbered
# 
# usage:
#    ./GenerateSubmitTarb.sh mypaper.tex myrefs.bib mypaper.bbl myrefs_loc.bib
#
# notes:
#     - fig names in the tex must have their extension, 
#       otherwise the variable "ext" has to be set by hand
#       all the figs are assumed with the same "ext"
#     - comments in the tex must be %%<space>text, 
#       otherwise change by hand "patCom"
# 
# sbernuz 20110818

# additional notes
#      bib file not necessary on arxiv
#      do not use patCom="%" because it deletes the \% 
#      

# store inputs into vars
texf=$1 # tex file
bibf=$2 # bib file 
bblf=$3 # bbl file
bibf_loc=$4 # bib file local

# output filenames
today=$(date +%Y%m%d)
newtexf=paper$today.tex
newbibf=refs$today.bib
newbblf=paper$today.bbl
newbibf_loc=refs_loc$today.bib

# patterns 
patBFig="\begin{figure"
patIFig="\includegraphics["
patEFig="\end{figure"
patCom="%% "

# create submission copies
cp $texf $newtexf
cp $bibf $newbibf
cp $bblf $newbblf
cp $bibf_loc $newbibf_loc

# extract tex fig and graphics lines
figl=$(fgrep -e "${patFig}" -e "${patGra}" $newtexf)

# counters
f=0
c=0
ab=(a b c d e f g h i l m n o p q r s t u v w x y z aa bb cc dd ee ff gg hh)

# for each fig ...
for l in $figl ; do 

    ##echo $l

    if [[ $l =~ "${patBFig}" ]]; then		

        # reset include fig counter
	c=0	

	# increment fig counter
	let f++       

	echo "==> begin fig no $f"

	continue

    fi

    if [[ $l =~ "${patIFig}" ]]; then

        # extract fig name and ext    
        #      \includegraphics[...]{$fname}
	fname=${l##*"]{"}
	fname=${fname%"}"*}	

	# rm path - no!
	##fname=${fname##*/}	

	# extract ext
	ext=${fname##*"."}
	##ext=".eps" #".pdf" # set by hand

        # generate new fig (cp)
	newfname=fig$(printf "%02d" $f)${ab[$c]} 	
	cp $fname $newfname.$ext	

        # replace tex fig-line with newfig name
        ##sed -ie 's%$fname%$newfname%g' $newtexf
        # older versions of sed:
	##sed 's#$fname#$newfname#g' $newtexf > tmp && mv tmp $newtexf
	sed "s%$fname%$newfname.$ext%g" $newtexf > tmp && mv tmp $newtexf

	echo "==* included fig no $f / $c :: $fname -> $newfname.$ext"
	
	let c++

	continue

    fi

    if [[ $l =~ "${patEFig}" ]]; then		

	if [ "$c" -eq 1 ]; then

    	    # correct previous deleting ab suf
	    oldfname=$newfname	    
	    newfname=fig$(printf "%02d" $f)

	    # generate new fig (mv)
	    mv $oldfname.$ext $newfname.$ext            
            
            # replace tex fig-line with newfig name
            ##sed -ie 's%$fname%$newfname%g' $newtexf
            # older versions of sed:
	    ##sed 's%\<$oldfname\>%$newfname%g' $newtexf > tmp && mv tmp $newtexf
    	    sed "s%$oldfname.$ext%$newfname.$ext%g" $newtexf > tmp && mv tmp $newtexf

	    echo "==* (single pic) fig no $f :: $oldfname.$ext -> $newfname.$ext"  
	    
	fi

	echo "<== end fig no $f"

	continue

    fi
    
done

##exit

# replace bib line
##sed -ie "s%${bibf%.bib}%${newbibf%.bib}%" $newtexf
# older versions of sed:
#sed "s%${bibf%.bib}%${newbibf%.bib}%g" $newtexf > tmp && mv tmp $newtexf
sed "s%${bibf%.bib,bibf_loc.bib}%${newbibf%.bib,newbibf_loc.bib}%g" $newtexf > tmp && mv tmp $newtexf

# remove partial line starting with comment pattern
sed "s/${patCom}.*//" $newtexf > tmp && mv tmp $newtexf

# generate tarball 
mkdir -vp submission$today
mv $newtexf $newbblf $newbibf $newbibf_loc fig*.$ext submission$today/ 
tar zcvf submission$today.tgz submission$today 

