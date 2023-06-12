NTER_CMAP = {
    "Z-type":"#7ad5a3",
    "X-type":"#322f26",
    "CoBaHMA-type":"#0babc1",
    "Y-type":"#b980d1",
    "Unknown-type":"#717171",
    "NA":"#717171"
}




// DOCUMENT IS READY 

$( document ).ready(function() {            
    $(window).trigger('resize');

    const strains = $.map(DATAS, function(v, i){            
        return i;
    });    

    SEQUENCES = {}
    $.each(DATAS,function(strainid,genomes){
        $.each(genomes,function(gid,datas){
            $.extend(SEQUENCES, datas["sequences"]);
        })
    })

    NTERS = []
    FLAGS = []
    $.each(SEQUENCES,function(seqid,seq){
        if ($.inArray( seq.nter , NTERS) == -1){
            NTERS.push(seq.nter)            
            $("#nter-select").append(
                `<option value='${seq.nter}'>${seq.nter}</option>`    
            )
        }
        if ($.inArray( seq.flag , FLAGS) == -1){
            FLAGS.push(seq.flag)            
            $("#flag-select").append(
                `<option value='${seq.flag}'>${seq.flag}</option>`    
            )
        }
    })

    generate_card()
    $("#overlay").hide()



    // FILTER CARD LIST WHEN ENTER KEY IS PRESSED
    $(".input_list").keypress(function(event) {
        if(event.which == 13){ // enter key
            var val = $(this).val().split(" ");
            var target = $(this).attr("ul")            
            $("#"+target + " li").show().filter(function() { 
                var li = $(this);
                var r = false;
                $.each(val,function(i,e){    
                                                 
                    if (li.children("a").attr("info").toLowerCase().trim().indexOf(e.toLowerCase()) == -1){                
                        r = true;
                     }             
                });                                
                return r;
            }).hide();
                                   
            $(".search-result").html(`<p style='font-style: italic;'>${$("#"+target + " li:visible").length} entries selected.</p>`)
        }  
        
      });

      $(".clear-input").click(function(){
        $(".input_list").val("")
        var target = $(this).attr("ul")        
        $("#"+target + " li").show()
        $(".search-result").html(`<p style='font-style: italic;'>${$("#"+target + " li:visible").length} entries selected.</p>`)
  })


      // ADD AND REMOVE FROM CART
      $("body").on("click",".add_to_cart",function(){
        var prot = $(this).attr("id") //.prev().attr("protein")
        var prot_in_cart = []
        $.each($("#Cart").children(),function(){
            prot_in_cart.push($(this).text())            
        })
        if ($.inArray( prot , prot_in_cart) == -1){
            $("#Cart").append(`<div class='cart_item'>${prot}</div>`)
        }
        $(`.add_to_cart[id="${prot}"]`).toggleClass("remove_from_cart")
    })

    $("body").on("click",".remove_from_cart",function(){
        var prot = $(this).attr("id") //.prev().attr("protein") 
        $.each($("#Cart").children(),function(){
            if ($(this).text()==prot ){
                $(this).remove()
            }
        })
    })

    $(".searchable_list").on("click",'.entry-card',function(){
        var seqid = $(this).attr("protein");
        var strain = $(this).attr("strain");
        var gid = $(this).attr("gid");
        render_entry(seqid,strain,gid);        
    })


    $(".sms_button").click(function(){
        target = $(this).attr("target")            
        
        $(".sms_button").removeClass("button_active")
        $(".sms_plot").hide()
        
        $(this).addClass("button_active")
        $("#"+target).fadeIn(10)
    })


    $(".oms_button").click(function(){
        target = $(this).attr("target")
        var id = $("#"+target).children().children().attr('id')//$(this).attr('id')
        
        var myPlot = document.getElementById(id);
        
        if (myPlot){
            myPlot.removeAllListeners('plotly_click')
            console.log("removing listeners")
        }
        
        $(".oms_button").removeClass("button_active")
        $(".oms_plot").hide()
        
        $(this).addClass("button_active")
        $("#"+target).fadeIn(10)
  
        if (myPlot && target!="sunburst" && target!="treemap"){
            
            myPlot.on('plotly_click', function(data){
        
                var entry_data = data.points.map(function(d){
                    return (d.data.customdata);
                });
                
                var seqid = entry_data[0][0].seqid;
                var strain = entry_data[0][1].Organism_Name;
                var gid = entry_data[0][2].Assembly;
                render_entry(seqid,strain,gid);           
            })
        }
       
    })

    $(".toggle-section").click(function(){
        $(this).toggleClass("toggle-off-section")
        if ($(this).hasClass("toggle-off-section")){
            $(this).parent().next().slideUp(1000)
        } else {
            $(this).parent().next().slideDown(1000)
        }

        
    })


    $("#Cards-container").on("click",".close-parent",function(){
        // $("#Cards").parent().slideUp(1000)       
        $("#Cards-container").slideUp(1000)       
    });


    $("#Cards").on("click",".fasta",function(){
        var seqid = $(this).attr("seq")
        var gid = $(this).attr("gid")
        var strain = $(this).attr("strain")
        var seqtype = $(this).attr("seqtype")
        var seq= DATAS[strain][gid]["sequences"][seqid]["sequence"]
        if (seqtype=="fna") {
            seq= DATAS[strain][gid]["sequences"][seqid]["ccyA_seq"]
        }
        var fasta = format_fasta(seqid, seq , DATAS[strain][gid]["sequences"][seqid])
        navigator.clipboard.writeText(fasta);
    
        var e = $(this);
        e.addClass("copied")
        setTimeout(function() {
            e.removeClass("copied")
        }, 2000);         
    })

    $("#Cards").on("click",".feature",function(){
        var seqid = $(this).attr("seq")
        var gid = $(this).attr("gid")
        var strain = $(this).attr("strain")
        var fid = $(this).attr("feature")
        console.log(fid)
        var feature= DATAS[strain][gid]["sequences"][seqid]["features"][fid]
        var identifier = seqid + "_" + feature.feature_id
        var seq = feature.feature_seq
        var fasta = format_fasta(identifier, seq , DATAS[strain][gid]["sequences"][seqid])
        navigator.clipboard.writeText(fasta);
    
        var e = $(this);
        e.addClass("copied")
        setTimeout(function() {
            e.removeClass("copied")
        }, 2000);         
    })



    $(".clear_cart").click(function(){
        $("#Cart").html("");
        $(".remove_from_cart").removeClass("remove_from_cart")
    });
});

//

function render_entry(seqid,strain,gid){
    console.log(seqid,strain,gid)
    render_genome(strain,gid);
    if (seqid){
        render_protein(strain,gid,seqid);
        render_gene(strain,gid,seqid);
        render_feature(strain,gid,seqid);
        render_hit(strain,gid,seqid);
    }

    $('html,body').animate({
        scrollTop: $("#Cards").offset().top},
        'slow');
}

// FUNCTION AND EVENT
function generate_card(){
    let div = "#card-ul"
    $.each( DATAS, function( strainid, genomes ){
            $.each( genomes ,function( gid, data ){
                var warning = "";
                if (Object.keys(data.sequences ).length >= 1 ){
                    
                    if (Object.keys(data.sequences ).length > 1){
              
                        warning =  "MULTIPLE";
                    }
                    
                    $.each(data.sequences , function(seqid,seqdatas){
                        
                        let infos = {"genotype":"ccyA-","flag":null,"nter":null,"cter":null,"warning":null}
                        
                        if (seqdatas.flag=="Calcyanin with known N-ter"){
                            infos["genotype"]="ccyA+";
                        } else {
                            infos["genotype"]="ccyA~";
                        }
                     
                        // let color="#ddd";
                        color = NTER_CMAP[seqdatas.nter];
                        infos["nter"] = seqdatas.nter;
                        infos["cter"] = seqdatas.cter;
                        infos["flag"] = seqdatas.flag;
                        infos["warning"] = warning;
    
                        let info_str = `${strainid} ${gid} ${infos.genotype} ${infos.nter} ${infos.cter} ${infos.flag} ${seqid} ${infos.warning}`
                        let e=`<li class='list-border'><a style='border-left:15px solid ${color};' 
                            class=entry-card info='${info_str}' protein='${seqid}' gid='${gid}' strain='${strainid}' >
                                <div class=li-flex>
                                    <div class='li-flex-item li-flex-item-large'>${strainid}</div>
                                    <div class=li-flex-item>${gid}</div>
                                    <div class=li-flex-item>${infos.genotype}</div>
                                    <div class=li-flex-item>${infos.flag}</div>
                                    <div class=li-flex-item>${infos.nter}</div>
                                    <div class=li-flex-item>${infos.cter}</div>
                                    <div class=li-flex-item>${seqid}</div>
                                    <div class=li-flex-item>${infos.warning}</div>
                                </div> 
                            </a>
                            <span class='add_to_cart icon' id='${seqid}'></span>
                            </li>`;
                        $(div).append(e)
                    })
                } else {
  
                    // strain-gid without calcyanin
                    let infos = {
                        "genotype":"ccyA-",
                        "flag":"",
                        "nter":"",
                        "cter":"",
                        "warning":"",
                    }
                    
                  
                    let info_str = `${strainid} ${gid} ${infos.genotype}  ${infos.cter} ${infos.flag}`
                    let color="#ddd";
                    let e=`<li class='list-border'><a style='border-left:15px solid ${color};' 
                    class=entry-card info='${info_str}' gid='${gid}' strain='${strainid}'>
                        <div class=li-flex>
                            <div class='li-flex-item li-flex-item-large' >${strainid}</div>
                            <div class=li-flex-item>${gid}</div>
                            <div class=li-flex-item>${infos.genotype}</div>
                            <div class=li-flex-item>${infos.flag}</div>
                            <div class=li-flex-item>${infos.nter}</div>
                            <div class=li-flex-item>${infos.cter}</div>
                            <div class=li-flex-item></div>
                            <div class=li-flex-item>${infos.warning}</div>
                        </div> 
                    </a>
                    </li>`;
                $(div).append(e)
                }
            })
        })
    }



function get_gid(strain) {
    return $.map(_DATAS_["strains"][strain], function(v, i){            
        return i;
    });            
}

function get_strain_from_acc(acc) {
    var strain_obj = null;
    $.each(_DATAS_.strains,function(straind, d){                
        $.each(d , function (dacc , _){                    
            if (dacc==acc){  
                strain_obj = _;                      
                return _;
            }
        })
    })            
    return strain_obj;
}



// RENDER GENOME AND PROTEIN    
function render_genome(strain,gid){
    $("#Cards-container").hide()
    $("#Genome-card").html("")        
    // $("#Protein-card").html("")
    
    d = DATAS[strain][gid]

    var url = "https://www.ncbi.nlm.nih.gov/data-hub/genome/" + gid
    
    content = `
    <div class='genome_card'>
        <div class=card_summary>
            <h1> Assembly </h1>
            <p><b>${strain}</b></p>
            <p>${gid}</p>
            <div class='button small_button inline-block'><a href='${url}' target='_blank' >NCBI</a></div>                    
            <div><span class="icon close-parent"></span></div>
        </div>
        <div class='key-val wrap'>
    ` // DIV NOT CLOSE ! 
    
    $.each(d,function(index,value){        
        if (index != "sequences" ){
            index = index[0].toUpperCase() + index.slice(1);
            content += render_shadow_key_val(index,value)
            // content+=`<div class=key-val-item><p class=key>${index}</p><p class=value>${value}</p></div>`
        }
        
    })    
    content = content + "</div></div>" // DIV CLOSE !
    $("#Genome-card").html(
        content
    )
    $("#Cards-container").slideDown(500)
}


function render_gene(strain,gid,seqid){
    seqrecord = DATAS[strain][gid]["sequences"][seqid];
  
    // $("#Protein-card").hide()
    if (seqid){
        content = `
        <div class='gene_card'>
            <div class=card_summary>
                <h1> Gene(s) </h1>
                <p><b>${seqid}</b></p>                                           
                <div></div>
            </div>
            <div class='key-val box-shadow'>
        ` // DIV NOT CLOSE !         
        content += render_key_val("Genomic region",seqrecord.ccyA_genomic_region)
        content += render_key_val("Start position",seqrecord.ccyA_start)
        content += render_key_val("End position",seqrecord.ccyA_stop)
        content += render_key_val("Frame",seqrecord.ccyA_frame)
        content += render_key_val("Partial",seqrecord.ccyA_partial)
        content += render_key_val("Pseudo",seqrecord.ccyA_pseudo)
        content += `
            <div class='copy fasta tooltip' seqtype=fna seq='${seqid}' strain='${strain}' gid='${gid}' >
                <span class="tooltiptext tooltiptext-left">
                    Copy ccyA gene in fasta format.
                </span>
            </div>`

        content = content + "</div></div>" // DIV CLOSE !
        $("#Gene-card").html(
            content
        )
    }
}


// 


    






function render_protein(strain,gid,seqid){
    seqrecord = DATAS[strain][gid]["sequences"][seqid];
  
    // $("#Protein-card").hide()
    if (seqid){
        let current_class = ""
        if ($(`li span.add_to_cart[id="${seqid}"]`).hasClass("remove_from_cart")){
            current_class = "remove_from_cart"
        }
        content = `
        <div class='protein_card'>
            <div class=card_summary>
                <h1> Protein(s) </h1>
                <p><b>${seqid}</b></p>

                <span class='add_to_cart icon ${current_class}' id='${seqid}'></span>
            </div>
            <div class='key-val box-shadow'>
        ` // DIV NOT CLOSE !         
        content += render_key_val("Calcyanin flag",seqrecord.flag)
        content += render_key_val("C-ter composition",seqrecord.cter)
        content += render_key_val("N-ter type",seqrecord.nter)
        content += render_key_val("Sequence length",seqrecord.sequence.length)
        content += render_key_val("Nearest neighbor",seqrecord.nter_neighbor)   
        content += `
        <div class='fasta copy tooltip' seqtype=faa seq='${seqid}' strain='${strain}' gid='${gid}' >
            <span class="tooltiptext tooltiptext-left">
                Copy calcyanin protein in fasta format.
            </span>
        </div>
        `     
        content = content + "</div></div>" // DIV CLOSE !
        $("#Protein-card").html(
            content
        )
    }
}



function render_feature(strain,gid,seqid){
    seqrecord = DATAS[strain][gid]["sequences"][seqid];
  
    // $("#Protein-card").hide()
    if (seqid){
        content = `
        <div class='feature_card'>
            <div class=card_summary>
                <h1> Feature(s) </h1>
                <p><b>${seqid}</b></p>                    
            </div>            
        ` // DIV NOT CLOSE !      

        $.each(seqrecord["features"],function(index,feature){      
            if (index % 2 == 0){
                color_class = "odd-color";
            } else {
                color_class = "even-color";
            }
            content += "<div class='key-val box-shadow "+color_class+"'>";
            content += _render_feature(index,feature);    
            content += `
                <div class='feature copy tooltip' feature='${index}'  seq='${seqid}' strain='${strain}' gid='${gid}'>
                    <span class="tooltiptext tooltiptext-left">
                        <p style='text-align:left;'>Copy ${feature.feature_id} for ${seqid} in fasta format.</p>
                    </span>
                </div>
                `
            content += "</div>";
        })

        content+="</div>"
        $("#Feature-card").html(
            content  
        )        
    }
}

function render_hit(strain,gid,seqid){
    seqrecord = DATAS[strain][gid]["sequences"][seqid];
    if (seqid){
        content = `
        <div class='hit_card'>
            <div class=card_summary>
                <h1> Raw Hits(s) </h1>
                <p><b>${seqid}</b></p>                    
            </div>            
        ` // DIV NOT CLOSE !      
        $.each(seqrecord["hits"],function(index,hit){      
            if (index % 2 == 0){
                color_class = "odd-color";
            } else {
                color_class = "even-color";
            }
            content += "<div class='key-val box-shadow "+color_class+"'>";
            content += _render_hit(index,hit);            
            content += "</div>";
        })

        content+="</div>"
        $("#Hit-card").html(
            content  
        )        
    }
}

function render_shadow_key_val(k,v){
    return "<div class='key-val-item box-shadow'><p class=key>"+k+"</p><p class=value>"+v+"</p></div>"
}

function render_key_val(k,v){
    return "<div class=key-val-item><p class=key>"+k+"</p><p class=value>"+v+"</p></div>"
}

function _render_feature(index,feature){
    feature_html = ""
    feature_html += render_key_val("Feature ID" , "<b>"+feature.feature_id+"</b>")
    feature_html += render_key_val("Feature src" , feature.feature_src)
    feature_html += render_key_val("Feature start" , feature.feature_start)
    feature_html += render_key_val("Feature end" , feature.feature_end)
    feature_html += render_key_val("Feature E-value" , feature["e-value"].toExponential())
    feature_html += render_key_val("Feature coverage" , feature.coverage)
    return feature_html
}

function _render_hit(index,hit){
    hit_html = ""
    hit_html += render_key_val("Hit ID" , "<b>"+hit.hit_src+"</b>")
    hit_html += render_key_val("Hit start" , hit.hit_start)
    hit_html += render_key_val("Hit end" , hit.hit_stop)
    hit_html += render_key_val("Hit E-value" ,hit.hit_e_value.toExponential())
    hit_html += render_key_val("Hit Coverage" , hit.hit_coverage)
    hit_html += render_key_val("Hit identity" , hit.hit_pident)
    return hit_html
}


function format_fasta(id,seq,seqO){
    var header = ">" + id + " " + seqO["nter"] + " " + seqO["sequence_src"] + " " + seqO["cter"]
    return header + "\n" + seq
}


function get_sequence(filter_key,filter_val){
    var sub = {}
    var sequences = {}
    $.each(_DATAS_["calcyanins"],function(seqid,seq_DATAS_){
        if (filter_val){                
            if (seq_DATAS_[filter_key] == filter_val ){
                sub[seqid]=seq_DATAS_
            } 
        } else {                
            sub[seqid]=seq_DATAS_                
        }
    })
    return sub
}

function get_feature_sequence(sequences,filter_key,filter_val){
    var seqs = {}

    $.each(sequences,function(seqid,seq_DATAS_){
        $.each(seq_DATAS_["features"],function(_,f_DATAS_){
            if (f_DATAS_[filter_key] == filter_val){
                seqs[seqid] = [filter_val, f_DATAS_.feature_seq];
            }
        })
    })
    return seqs
}    
    

function downloadTxt(txt,filename) {
    var blob = new Blob([txt], {
        type: "text/plain;charset=utf-8"
    });
    var a = document.createElement("a");
    a.href = window.URL.createObjectURL(blob, {type: "text/plain"});
    a.download = filename;
    a.click();
};
    
   
$(".expand_next").click(function(){
    $(this).parent().next().toggleClass("constrain_height_div")
})

function get_feature_seq(seqid, seqdatas,feature_id){
    let sequences = []
    let cpt = 0;
    $.each(seqdatas.features, function(fid , feature){    
        if (fid == feature_id){
            cpt += 1;
            let fasta = `>${seqid} ${feature_id} feature_no_${cpt}\n${feature.feature_seq}`;
            sequences.push(fasta)                
        }
    })
    return sequences
}


function get_seq_in_fasta_format(seqid,ftype){
    let seq_datas = SEQUENCES[seqid]
    let header = `>${seqid} ${seq_datas.nter} ${seq_datas.cter}  ${seq_datas.flag} ${seq_datas.sequence_src}` 
    let fasta = [];
    if (ftype=="faa"){
        fasta = [`${header}\n${seq_datas.sequence}`]
    } else if (ftype == "fna"){
        fasta = [`${header}\n${seq_datas.ccyA_seq}`]
    } else {
        fasta = get_feature_seq(seqid,seq_datas,ftype)
    }
    return fasta    
}

$(".download-cart").click(function(){
    let ftype = $(this).prev().find(":selected").val();
    let sequences_in_cart = []
    if (ftype){
        $("#Cart").children().each(function(){
            sequences_in_cart.push($(this).text())
        })
        if (sequences_in_cart.length>0){
            let fasta = [];
            $.each(sequences_in_cart,function(_,seqid){
                fasta = $.merge(fasta,get_seq_in_fasta_format(seqid,ftype));        
            })
            downloadTxt(fasta.join("\n"),`cart_${ftype}.fasta`)
        } else {
            $(this).next().css("display","inline-block")
            $(this).next().html("<p style='font-style: italic;'>Cart is empty.<p>")
            $(this).next().delay(3000).fadeOut('slow');
        }
    }
})

$(".download-batch").click(function(){
    let nter_subset = $(this).prev().prev().prev().find(":selected").val();
    let flag_subset = $(this).prev().prev().find(":selected").val();
    let ftype = $(this).prev().find(":selected").val();
    let subset = [];
    // first filter : based on their n-ter value
    subset = $.map(SEQUENCES, function(v, i){   
        if (nter_subset=="all"){
            return i;
        } else {
            if (SEQUENCES[i].nter==nter_subset){
                return i;
            }    
        }           
    });
    // second filter :  based on their flag value    
    subset = $.map(subset, function(i, _){           
        if (flag_subset=="all"){
            return i;
        } else {
            if (SEQUENCES[i].flag==flag_subset){
                return i;
            }    
        }           
    });

    if (ftype){
        if (subset.length>0){
            let fasta = [];
            $.each(subset,function(_,seqid){
                fasta = $.merge(fasta,get_seq_in_fasta_format(seqid,ftype));        
            })
           
            downloadTxt(fasta.join("\n"),`batch_${ftype}.fasta`)
        } else {
            $(this).next().css("display","inline-block")
            $(this).next().html("<p style='font-style: italic;'>Nothing match your filter.</p>")
            $(this).next().delay(3000).fadeOut('slow');
        }
    }
})

