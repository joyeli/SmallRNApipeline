const d3 = require( "d3" );

module.exports = {

    drewMiRNA: function( miRNA, label, left ){
        const mirnaLabel = label.selectAll( ".miRNA-label-" + miRNA )
            .data([ 0 ])
            .enter().append( "g" )
            .attr( "class", "miRNA-label-" + miRNA )
            .attr( "font-family", "sans-serif" )
            .attr( "font-size", 12 )
            .attr( "transform", "translate(-" + left + ",-14)" )
            .append( "text" )
            .attr( "text-anchor", "start" )
            .text( miRNA );
    },

    drewLabelx: function( labelsSeed, file_name, arm, label, height, x, y ) {
        const columnLabels = label.selectAll( ".column-label-" + file_name + arm )
            .data( labelsSeed )
            .enter().append( "g" )
            .attr( "class", "column-label-" + file_name + arm )
            .attr( "font-family", "sans-serif" )
            .attr( "font-size", 12 )
            .attr( "transform", function( d, i ){ return "translate(" + x(i) + "," + height + ")"; });

        columnLabels.append( "line" )
            .style( "stroke", "black" )
            .style( "stroke-width", "1px" )
            .attr(  "x1", x.bandwidth() / 2 )
            .attr(  "x2", x.bandwidth() / 2 )
            .attr(  "y1", 0 )
            .attr(  "y2", 5 );

        columnLabels.append( "text" )
            .attr( "x", 0 )
            .attr( "y", y.bandwidth() / 2 )
            .attr( "dx", ".6em" )
            .attr( "dy", ".9em" )
            .attr( "text-anchor", "end" )
            .attr( "transform", "rotate(-70)" )
            .text( function( d, i ){ return d; });

        const headerLabel = label.selectAll( ".header-label-" + file_name + arm )
            .data([ 0 ])
            .enter().append( "g" )
            .attr( "class", "header-label-" + file_name + arm)
            .attr( "font-family", "sans-serif" )
            .attr( "font-size", 12 )
            .attr( "transform", "translate(" + (x(0)+(x.step()*(labelsSeed.length/2))) + ",-4)" )
            .append( "text" )
            .attr( "text-anchor", "middle" )
            .text( arm );
    },

    drewHeader: function( labelsSeed, file_name, label, argv, x, y ) {
        const headerLabel = label.selectAll( ".header-label-" + file_name )
            .data([ 0 ])
            .enter().append( "g" )
            .attr( "class", "header-label-" + file_name )
            .attr( "font-family", "sans-serif" )
            .attr( "font-size", 12 )
            .attr( "transform", "translate(" + ( x(0) + ( x.step() * ((
                                        argv.arm == "5p" ? labelsSeed[ "5p" ].length :
                                        argv.arm == "3p" ? labelsSeed[ "3p" ].length :
                                        labelsSeed[ "5p" ].length + labelsSeed[ "3p" ].length ) /2 ))) + ",-14 )" )
            .append( "text" )
            .attr( "text-anchor", "middle" )
            .text( file_name );
    },

    drewLabely: function( labelsLength, label, y ) {
        const rowLabels = label.selectAll( ".row-label" )
            .data( labelsLength )
            .enter().append( "g" )
            .attr( "class", "row-label" )
            .attr( "font-family", "sans-serif" )
            .attr( "font-size", 12 )
            .attr( "transform", function( d, i ){ return "translate(" + 0 + "," + y(i) + ")"; });

        rowLabels.append( "line" )
            .style( "stroke", "black" )
            .style( "stroke-width", "1px" )
            .attr(  "x1", 0 )
            .attr(  "x2", -5 )
            .attr(  "y1", y.bandwidth() / 2 )
            .attr(  "y2", y.bandwidth() / 2 );

        rowLabels.append( "text" )
            .attr( "x", -8 )
            .attr( "y", y.bandwidth() / 2 )
            .attr( "dy", ".32em" )
            .attr( "text-anchor", "end" )
            .text( function( d, i ){ return d; });
    },

    drewScalebar: function( type, svg, height, legendMargin, widthLegend, startColor, endColor, minValue, maxValue ){
        const key = svg.append( "g" ).attr( "transform", "translate(" + legendMargin + ",0)" );

        const legend = key
            .append( "defs" )
            .append( "linearGradient" )
            .attr(   "id", "gradient" + type )
            .attr(   "x1", "100%" )
            .attr(   "y1", "0%" )
            .attr(   "x2", "100%" )
            .attr(   "y2", "100%" )
            .attr(   "spreadMethod", "pad" );

        legend.append( "stop" )
            .attr( "offset", "0%" )
            .attr( "stop-color", endColor )
            .attr( "stop-opacity", 1 );

        legend.append( "stop" )
            .attr( "offset", "100%" )
            .attr( "stop-color", startColor )
            .attr( "stop-opacity", 1 );

        key.append( "rect" )
            .attr(  "width", widthLegend/2 )
            .attr(  "height", height )
            .style( "fill", "url(#gradient" + type + ")" )
            .attr(  "transform", "translate(0,0)" );

        const scale_y = d3.scaleLinear().range([ height, 0 ]).domain([ minValue, maxValue ]);
        const yAxis = d3.axisRight().scale( scale_y );

        key.append( "g" )
            .attr(  "class", "y axis" )
            .attr(  "transform", "translate(41,0)" )
            .call(  yAxis );

        const gap = svg.selectAll( ".gap" )
            .data([ 0 ])
            .enter().append( "line" )
            .attr( "class", ".gap" )
            .attr( "transform", "translate(" + (legendMargin-5) + ",-15)" )
            .style( "stroke", "slategray" )
            .style( "stroke-width", "1px" )
            .attr( "x1", 0 ) 
            .attr( "x2", 0 ) 
            .attr( "y1", 0 ) 
            .attr( "y2", height + 80 ); 
    }
};
