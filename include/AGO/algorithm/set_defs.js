const jsdom = require( "jsdom" );
const d3 = require( "d3" );

const seedFuns = require( "./seed_funs" );
const dataParser = require( "./data_parser" );

module.exports = {

    get_defs: function( argv ) {

        const { JSDOM } = jsdom;
        const dom = new JSDOM();
        const document = dom.window.document;

        const files      = dataParser.readFiles( argv );
        const seed_index = seedFuns.getSeedIndex( files, argv );
        const labelsSeed = seedFuns.getSeedLabel( seed_index, argv );

        let datasize = 0;

        for( let file_name in files ){ datasize++; }

        const margin = {
            top    : 25,
            bottom : 60,
            left   : 40,
            right  : 10,
            gap    : 10
        };

        const height      = 615;
        const widthCell   = 36;
        const widthLegend = 90;

        const chartHeight = height + margin.top + margin.bottom;
        const chartWidth  = ( widthCell * datasize * (
                        argv.arm == "5p" ? labelsSeed[ "5p" ].length :
                        argv.arm == "3p" ? labelsSeed[ "3p" ].length :
                        labelsSeed[ "5p" ].length + labelsSeed[ "3p" ].length ))
                + ( margin.gap * datasize ) + margin.right +
                ( widthLegend * ( argv.mode == "heatmap" ? 2 : 1 ));

        const labelsLength = [ "Total",
            "30", "29", "28", "27", "26", "25", "24", "23",
            "22", "21", "20", "19", "18", "17", "16", "15"
                ];

        const numRows = labelsLength.length;

        const svg = d3.select( document.body ).append( "svg" )
            .attr( "xmlns", "http://www.w3.org/2000/svg" )
            .attr( "version", "1.1" )
            .attr( "width",  chartWidth )
            .attr( "height", chartHeight )
            .append( "g" )
            .attr( "transform", "translate(" + margin.left + "," + margin.top + ")" );

        const label = svg.append( "g" ).attr( "class", "labels" );

        const y = d3.scaleBand().domain( d3.range( numRows )).rangeRound([ 0, height ]);

        return {

            document     : document,
            svg          : svg,
            label        : label,
            y            : y,

            files        : files,
            seed_index   : seed_index,
            labelsSeed   : labelsSeed,

            labelsLength : labelsLength,
            datasize     : datasize,

            startColor   : { value: "white",      density: "white"   },
            endColor     : { value: "dodgerblue", density: "dimgrey" },

            numCols      : { "5p": labelsSeed[ "5p" ].length, "3p": labelsSeed[ "3p" ].length },
            numRows      : numRows,

            margin       : margin,
            height       : height,

            widthCell    : widthCell,
            widthLegend  : widthLegend,

            chartHeight  : chartHeight,
            chartWidth   : chartWidth,

            legendMargin : margin.gap + ( widthCell * datasize * (
                        argv.arm == "5p" ? labelsSeed[ "5p" ].length :
                        argv.arm == "3p" ? labelsSeed[ "3p" ].length :
                        labelsSeed[ "5p" ].length + labelsSeed[ "3p" ].length ))
        }
    }
};
