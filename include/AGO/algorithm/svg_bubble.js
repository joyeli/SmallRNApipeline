const d3 = require( "d3" );
module.exports = {

    tails2path: function( seed_index, tail_data, rwth, rdis, islog2 ){

        let tail_path = [];
        let empty_array = [];

        const empty_object = {
                GMPM       : 0,
                GM         : 0,
                PM         : 0,
                A_Tail     : 0,
                C_Tail     : 0,
                G_Tail     : 0,
                T_Tail     : 0,
                Other_Tail : 0
        };

        for( let seed in seed_index ){
            empty_array.push( empty_object );
        }

        empty_array.push( empty_object );
        tail_path.push( empty_array );

        for( let idx = 15; idx >= 0; idx-- ){

            let row_array = [];
            for( let seed in seed_index ){

                if( !( seed in tail_data )){
                    row_array.push( empty_object )
                } else {
                    row_array.push({
                        GMPM       : rdis *( Math.abs( islog2 ? Math.log2( tail_data[ seed ][ "GMPM" ][ idx ]) : tail_data[ seed ][ "GMPM" ][ idx ])),
                        GM         : rdis *( Math.abs( islog2 ? tail_data[ seed ][ "GM" ][ idx ] / tail_data[ seed ][ "GMPM" ][ idx ] * Math.log2( tail_data[ seed ][ "GMPM" ][ idx ]) : tail_data[ seed ][ "GM"   ][ idx ])),
                        PM         : rdis *( Math.abs( islog2 ? tail_data[ seed ][ "PM" ][ idx ] / tail_data[ seed ][ "GMPM" ][ idx ] * Math.log2( tail_data[ seed ][ "GMPM" ][ idx ]) : tail_data[ seed ][ "PM"   ][ idx ])),
                        A_Tail     : rwth *( tail_data[ seed ][ "A_Tail"     ][ idx ] / tail_data[ seed ][ "PM" ][ idx ] ),
                        C_Tail     : rwth *( tail_data[ seed ][ "C_Tail"     ][ idx ] / tail_data[ seed ][ "PM" ][ idx ] ),
                        G_Tail     : rwth *( tail_data[ seed ][ "G_Tail"     ][ idx ] / tail_data[ seed ][ "PM" ][ idx ] ),
                        T_Tail     : rwth *( tail_data[ seed ][ "T_Tail"     ][ idx ] / tail_data[ seed ][ "PM" ][ idx ] ),
                        Other_Tail : rwth *( tail_data[ seed ][ "Other_Tail" ][ idx ] / tail_data[ seed ][ "PM" ][ idx ] )
                    })
                };
            }

            row_array.push( empty_object );
            tail_path.push( row_array );
        }
        return tail_path;
    },

    drewBubble: function( data, file_name, type, arm, svg, height, x, y, seed_index, tail_data, maxValue, islog2 ) {

        let rdis = y.bandwidth() / maxValue / 2;
        let tail_path = module.exports.tails2path( seed_index, tail_data, y.bandwidth()/2, rdis, islog2 );

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
            .attr( "transform", function( d, i ){ return "translate(" + x(i) + ", 0)"; });

        let j = -1;
        let jend = tail_path.length -1;

        cell.append( "path" )
            .style( "fill", "red" )
            .attr( "transform", "translate(" + x.bandwidth()/2 + "," + y.bandwidth()/2 + ")" )
            .attr( "d", d3.arc()
                .startAngle( -.2 * Math.PI )
                .endAngle( .2 * Math.PI )
                .innerRadius( 0 )
                .outerRadius( function( d, i ){
                    if( i == 0 ){ j++; }
                    if( j > jend ){ j = 0; }
                    return tail_path[j][i] != null ? tail_path[j][i][ "A_Tail" ] : 0;
                }));

        cell.append( "path" )
            .style( "fill", "blue" )
            .attr( "transform", "translate(" + x.bandwidth()/2 + "," + y.bandwidth()/2 + ")" )
            .attr( "d", d3.arc()
                .startAngle( .2 * Math.PI )
                .endAngle( .6 * Math.PI )
                .innerRadius( 0 )
                .outerRadius( function( d, i ){
                    if( i == 0 ){ j++; }
                    if( j > jend ){ j = 0; }
                    return tail_path[j][i] != null ? tail_path[j][i][ "C_Tail" ] : 0;
                }));

        cell.append( "path" )
            .style( "fill", "darkorange" )
            .attr( "transform", "translate(" + x.bandwidth()/2 + "," + y.bandwidth()/2 + ")" )
            .attr( "d", d3.arc()
                .startAngle( .6 * Math.PI )
                .endAngle( 1 * Math.PI )
                .innerRadius( 0 )
                .outerRadius( function( d, i ){
                    if( i == 0 ){ j++; }
                    if( j > jend ){ j = 0; }
                    return tail_path[j][i] != null ? tail_path[j][i][ "G_Tail" ] : 0;
                }));

        cell.append( "path" )
            .style( "fill", "green" )
            .attr( "transform", "translate(" + x.bandwidth()/2 + "," + y.bandwidth()/2 + ")" )
            .attr( "d", d3.arc()
                .startAngle( 1 * Math.PI )
                .endAngle( 1.4 * Math.PI )
                .innerRadius( 0 )
                .outerRadius( function( d, i ){
                    if( i == 0 ){ j++; }
                    if( j > jend ){ j = 0; }
                    return tail_path[j][i] != null ? tail_path[j][i][ "T_Tail" ] : 0;
                }));

        cell.append( "path" )
            .style( "fill", "gray" )
            .attr( "transform", "translate(" + x.bandwidth()/2 + "," + y.bandwidth()/2 + ")" )
            .attr( "d", d3.arc()
                .startAngle( 1.4 * Math.PI )
                .endAngle( 1.8 * Math.PI )
                .innerRadius( 0 )
                .outerRadius( function( d, i ){
                    if( i == 0 ){ j++; }
                    if( j > jend ){ j = 0; }
                    return tail_path[j][i] != null ? tail_path[j][i][ "Other_Tail" ] : 0;
                }));

        cell.append( "path" )
            .style( "fill", "black" )
            .attr( "transform", "translate(" + x.bandwidth()/2 + "," + y.bandwidth()/2 + ")" )
            .attr( "d", d3.arc()
                .startAngle( 0 * Math.PI )
                .endAngle( 2 * Math.PI )
                .innerRadius( function( d, i ){
                    if( i == 0 ){ j++; }
                    if( j > jend ){ j = 0; }
                    return tail_path[j][i] != null ? tail_path[j][i][ "GMPM" ] : 0;
                })
                .outerRadius( function( d, i ){
                    return tail_path[j][i] != null ? tail_path[j][i][ "GMPM" ]+1 : 0;
                }));

        cell.append( "circle" )
            .attr( "transform", "translate(" + x.bandwidth()/2 + "," + y.bandwidth()/2 + ")" )
            .style( "stroke-dasharray", ( ".5, .5" ))
            .style( "stroke-width", 1 )
            .style( "stroke", "black" )
            .style( "fill", "none" )
            .attr(  "cx", 0 )
            .attr(  "cy", 0 )
            .attr(  "r", function( d, i ){
                if( i == 0 ){ j++; }
                if( j > jend ){ j = 0; }
                return tail_path[j][i] != null ? tail_path[j][i][ "PM" ]+1 : 0 });

        row.selectAll( ".cell-" + file_name + type + arm )
            .data( function( d, i ){ return data[i]; })
            .attr( "id", function( d ){ return d == "" || d == null ? "remove" : "cell"; });

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
