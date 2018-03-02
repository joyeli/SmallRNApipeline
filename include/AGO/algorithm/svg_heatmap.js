const d3 = require( "d3" );
module.exports = {

    drewHeatmap: function( data, file_name, type, arm, svg, height, x, y, startColor, endColor, minValue, maxValue ) {
        const colorMap = d3.scaleLinear().domain([ minValue, maxValue ]).range([ startColor, endColor ]);

        const row = svg.selectAll( ".row-" + file_name + type + arm )
            .data( data )
            .enter().append( "g" )
            .attr( "class", "row-" + file_name + type + arm )
            .attr( "transform", function( d, i ){ return "translate(0," + y(i) + ")"; })
            .attr( "id", function( d ){ let cell_count = 0;
                for( let j = 0, jlen = d.length; j < jlen; j++ ){
                    if( d[j] == "" ){ cell_count++; }
                }
                return cell_count == d.length ? "remove" : "row";
            });

        const cell = row.selectAll( ".cell-" + file_name + type + arm )
            .data( function( d ){ return d; })
            .enter().append( "g" )
            .attr( "class", "cell-" + file_name + type + arm )
            .attr( "font-family", "sans-serif" )
            .attr( "font-size", 8 )
            .attr( "transform", function( d, i ){ return "translate(" + x(i) + ", 0)"; });

        cell.append( "rect" )
            .attr(  "width", x.bandwidth() )
            .attr(  "height", y.bandwidth() )
            .style( "stroke-width", 0 );

        cell.append( "text" )
            .attr(  "dy", ".32em")
            .attr(  "x", x.bandwidth() / 2 )
            .attr(  "y", y.bandwidth() / 2 )
            .attr(  "text-anchor", "middle" )
            .style( "fill", function( d, i ){ return d >= maxValue/2 ? 'white' : 'black'; })
            .text( function( d, i ){ return d == "" || d == null ? "" : d.toPrecision(5); });

        row.selectAll( ".cell-" + file_name + type + arm )
            .data( function( d, i ){ return data[i]; })
            .attr( "id", function( d ){ return d == "" || d == null ? "remove" : "cell"; })
            .style( "fill", colorMap );

        svg.selectAll( "#remove" ).remove();

        if( type == "value" ){
            const gap = svg.selectAll( ".gap-" + file_name + arm )
                .data([ 0 ])
                .enter().append( "line" )
                .attr( "class", ".gap-" + file_name + arm )
                .attr( "transform", "translate(" + (x(0)-5) + ",-15)" )
                .style( "stroke", arm == "5p" ? "slategray" : "lightgray" )
                .style( "stroke-width", "1px" )
                .attr( "x1", 0 ) 
                .attr( "x2", 0 ) 
                .attr( "y1", 0 ) 
                .attr( "y2", arm == "5p" ? height + 80 : height + 50 ); 
        }
    }
};
