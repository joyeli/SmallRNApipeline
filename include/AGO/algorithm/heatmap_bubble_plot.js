const fs = require( "fs" );
const d3 = require( "d3" );

const setDefs  = require( "./set_defs" );
const dataParser = require( "./data_parser" );

const svgLabel = require( "./svg_label" );
const svgHeatmap = require( "./svg_heatmap" );
const svgBubble = require( "./svg_bubble" );

const makeTable = require( "./table_funs" );

const argv = require( "yargs" )
    .option({
        "input": {
            alias: "i",
            type : "string",
            describe: "Input miRNA name w/out seed, like MIR17HG-206-3p",
            demandOption: true
        },
        "rename": {
            alias: "r",
            type : "string",
            describe: "Rename miRNA name w/out arm and seed, like miR-92a"
        },
        "log2": {
            alias: "l",
            type : "boolean",
            default:  false,
            describe: "Set Log2 on"
        },
        "type": {
            alias: "t",
            choices: [ "GMPM", "GM", "PM", "A_Tail", "C_Tail", "G_Tail", "T_Tail", "Other_Tail" ],
            default: "GMPM",
            describe: "Select a type of expression data, Tails are only work with mode[heatmap]"
        },
        "arm": {
            alias: "a",
            choices: [ "5p3p", "5p", "3p" ],
            default: "5p3p",
            describe: "Select the arm to display"
        },
        "arm5seq": {
            alias: "5",
            type : "string",
            default: "",
            describe: "Input raw sequence of 5p reads"
        },
        "arm3seq": {
            alias: "3",
            type : "string",
            default: "",
            describe: "Input raw sequence of 3p reads"
        },
        "minlen": {
            alias: "n",
            type : "number",
            default: 15,
            describe: "Input min length for plot axis"
        },
        "maxlen": {
            alias: "x",
            type : "number",
            default: 30,
            describe: "Input max length for plot axis"
        },
        "mode": {
            alias: "m",
            choices: [ "heatmap", "bubble" ],
            default: "heatmap",
            describe: "Select mode of the chart"
        },
        "files": {
            alias: "f",
            type : "array",
            describe: "Input files from AGO data",
            demandOption: true
        }
    }).argv;

const defines = setDefs.get_defs( argv );

let datas = {};
let maxArray = { value: Array(), density: Array() };
let minArray = { value: Array(), density: Array() };

for( let file_name in defines.files ){
    datas[  file_name ] = dataParser.getData( argv, defines.files[ file_name ], defines.seed_index );

    for( let arm in { "5p":0, "3p":0 }){
        if( argv.arm != "5p3p" ){ if( arm != argv.arm ){ continue; }}

        for( let type in maxArray ){
            maxArray[ type ].push( d3.max( datas[ file_name ][ type ][ arm ], function(l){ return d3.max( l, function(d){ return d; }); }));
            minArray[ type ].push( d3.min( datas[ file_name ][ type ][ arm ], function(l){ return d3.min( l, function(d){ return d; }); }));
        }
    }
}

const maxValue = {
    value:   d3.max( maxArray[ "value"   ], function( d ){ return d; }),
    density: d3.max( maxArray[ "density" ], function( d ){ return d; })
};

const minValue = {
    value:   d3.min( minArray[ "value"   ], function( d ){ return d; }),
    density: d3.min( minArray[ "density" ], function( d ){ return d; })
};

let bandSize = 0;
let bandArray = [];

for( let file_name in datas ){
    for( let arm in datas[ file_name ][ "value" ] ){
        if( argv.arm != "5p3p" ){ if( arm != argv.arm ){ continue; }}

        bandSize += defines.labelsSeed[ arm ].length;
        bandArray.push( bandSize );
    }
}

let x = 0;
let count = 0;

for( let file_name in datas ){

    x = d3.scaleBand().domain( d3.range(
                argv.arm == "5p" ? defines.numCols[ "5p" ] :
                argv.arm == "3p" ? defines.numCols[ "3p" ] :
                defines.numCols[ "5p" ] + defines.numCols[ "3p" ] ))
        .rangeRound([
            count == 0 ? defines.margin.gap : defines.margin.gap + ( defines.widthCell * bandArray[ count -1 ]),
            defines.widthCell * bandArray[ argv.arm == "5p" ? count : argv.arm == "3p" ? count : count +1 ]
        ]);

    svgLabel.drewHeader( defines.labelsSeed, file_name, defines.label, argv, x, defines.y );

    for( let arm in datas[ file_name ][ "value" ] ){

        if( argv.arm != "5p3p" ){ if( arm != argv.arm ){ continue; }}

        x = d3.scaleBand().domain( d3.range( defines.numCols[ arm ] )).rangeRound([
                count == 0 ? defines.margin.gap : defines.margin.gap + ( defines.widthCell * bandArray[ count -1 ]),
                defines.widthCell * bandArray[ count ]
                ]);

        if( argv.mode == "heatmap" ){
            let type = "value";

            svgHeatmap.drewHeatmap( datas[ file_name ][ type ][ arm ], file_name, type, arm, defines.svg, defines.height, x, defines.y,
                    defines.startColor[ type ], defines.endColor[ type ], minValue[ type ], maxValue[ type ] );
        }

        if( argv.mode == "bubble" ){
            let type = "value";

            svgBubble.drewBubble( datas[ file_name ][ type ][ arm ], file_name, type, arm, defines.svg, defines.height, x, defines.y,
                    defines.seed_index[ arm ], defines.files[ file_name ][ arm ], maxValue[ type ], argv.log2 );
        }

        let type = "density";
        svgHeatmap.drewHeatmap( datas[ file_name ][ type ][ arm ], file_name, type, arm, defines.svg, defines.height, x, defines.y,
                defines.startColor[ type ], defines.endColor[ type ], minValue[ type ], maxValue[ type ] );

        svgLabel.drewLabelx( defines.labelsSeed[ arm ], file_name, arm, defines.label, defines.height, x, defines.y );
        count++;
    }
}

let scale_count = 0;

for( let type in maxValue ){
    if( argv.mode != "heatmap" && type == "value" ){ continue; }

    svgLabel.drewScalebar( type, defines.svg, defines.height, ( defines.legendMargin + defines.widthLegend * scale_count ), defines.widthLegend,
            defines.startColor[ type ], defines.endColor[ type ], minValue[ type ], maxValue[ type ] );

    scale_count++;
}

svgLabel.drewMiRNA( argv.rename == null ? argv.input : argv.rename, defines.label, defines.margin.left );
svgLabel.drewLabely( defines.labelsLength, defines.label, defines.y );

console.log( defines.document.body.innerHTML );
makeTable.log( defines.files, datas, argv );
