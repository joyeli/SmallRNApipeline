const fs = require( "fs" );
const split = require( "split-string" );

module.exports = {

    readFiles: function( argv ) {

        let files = {};
        let seed_obj = { "5p" : argv.arm5seq, "3p" : argv.arm3seq };

        for( let i = 0, ilen = argv.files.length; i < ilen; i++ ){
            files[ split( argv.files[i], "." )[0] ] = module.exports.readData( argv.files[i], argv.input, argv.minlen, argv.maxlen, seed_obj );
        }

        return files;
    },

    readData: function( filein, miRNA, minlen, maxlen, seed_obj ) {
        // miRNA_Length	A_Tail	C_Tail	G_Tail	U_Tail	Other_Tail	GM
        // MIR17HG-206-3p_ATTGCAC*:15	0	0	0	0	0	1.27017

        let idx = 0;
        let lens = {};

        for( let i = minlen; i <= maxlen; ++i ){
            lens[ i ] = idx;
            idx++;
        }

        let miRNA_5p = miRNA + "-5p";
        let miRNA_3p = miRNA + "-3p";

        let file = { "5p": {}, "3p":{} };
        let indata_split_line = split( fs.readFileSync( filein, "utf8" ), "\n" );
        let indata_split_head = split( indata_split_line[0], "\t" );

        for( let i = 1, ilen = indata_split_line.length -1; i < ilen; i++ ){
            let temp_miRNA = split( indata_split_line[i], "_" );

            if( temp_miRNA[0] != miRNA_5p && temp_miRNA[0] != miRNA_3p ){ continue; }

            let arm_split = split( temp_miRNA[0], "-" );
            let arm  = arm_split[ arm_split.length -1 ];
            let seed = String( split( temp_miRNA[1], ":" )).substr( 0, 7 );

            if( seed_obj[ arm ].indexOf( seed ) == -1 ) continue;

            let splits = split( indata_split_line[i], "\t" );
            let length = split( splits[0], ":" )[1];

            if( !( length in lens )) continue;

            if( file[ arm ][ seed ] == null ){
                file[ arm ][ seed ] = {
                    GM:         Array( lens.length ).fill( 0 ),
                    PM:         Array( lens.length ).fill( 0 ),
                    GMPM:       Array( lens.length ).fill( 0 ),
                    A_Tail:     Array( lens.length ).fill( 0 ),
                    C_Tail:     Array( lens.length ).fill( 0 ),
                    G_Tail:     Array( lens.length ).fill( 0 ),
                    U_Tail:     Array( lens.length ).fill( 0 ),
                    Other_Tail: Array( lens.length ).fill( 0 )
                };
            }

            for( let j = 1, jlen = indata_split_head.length; j < jlen; j++ ){
                file[ arm ][ seed ][ indata_split_head[j] ][ lens[ length ]] =
                (
                    file[ arm ][ seed ][ indata_split_head[j] ][ lens[ length ]] != null ?
                  ( file[ arm ][ seed ][ indata_split_head[j] ][ lens[ length ]] + Number( splits[j] ))
                  : Number( splits[j] )
                );
            }

            file[ arm ][ seed ][ "PM" ][ lens[ length ]]
                = file[ arm ][ seed ][ "A_Tail" ][ lens[ length ]]
                + file[ arm ][ seed ][ "C_Tail" ][ lens[ length ]]
                + file[ arm ][ seed ][ "G_Tail" ][ lens[ length ]]
                + file[ arm ][ seed ][ "U_Tail" ][ lens[ length ]]
                + file[ arm ][ seed ][ "Other_Tail" ][ lens[ length ]];

            file[ arm ][ seed ][ "GMPM" ][ lens[ length ]]
                = file[ arm ][ seed ][ "GM" ][ lens[ length ]]
                + file[ arm ][ seed ][ "PM" ][ lens[ length ]];
        }

        return file;
    },

    getData: function( argv, file, seed_index ){

        let lens = [];
        for( let i = argv.minlen; i <= argv.maxlen; ++i ){
            lens.push( i );
        }

        let file_data  = {
            "value"   : { "5p": Array( lens.length ), "3p": Array( lens.length ) },
            "density" : { "5p": Array( lens.length ), "3p": Array( lens.length ) }
        };

        for( let arm in seed_index ){
            if( argv.arm != "5p3p" ){ if( arm != argv.arm ){ continue; }}

            let total_array = Array();
            total_array[0] = 0;

            for( let length = 0; length < lens.length; length++ ){
                let values    = Array();
                let densities = Array();
                let total = 0;

                for( let seed in seed_index[ arm ] ){
                    values   [ seed_index[ arm ][ seed ]] = "";
                    densities[ seed_index[ arm ][ seed ]] = "";

                    total_array[ seed_index[ arm ][ seed ]] =
                        total_array[ seed_index[ arm ][ seed ]] != null
                        ? total_array[ seed_index[ arm ][ seed ]]
                        : 0;

                    if( file[ arm ][ seed ] != null ){

                        total_array[ seed_index[ arm ][ seed ]] += file[ arm ][ seed ][ argv.mode == "heatmap" ? argv.type : "GMPM" ][ length ];
                        total += file[ arm ][ seed ][ argv.mode == "heatmap" ? argv.type : "GMPM" ][ length ];

                        values[ seed_index[ arm ][ seed ] ] = file[ arm ][ seed ][ argv.mode == "heatmap" ? argv.type : "GMPM" ][ length ] == 0 ? "" :
                            argv.log2 ? Math.log2( file[ arm ][ seed ][ argv.mode == "heatmap" ? argv.type : "GMPM" ][ length ] ) :
                            file[ arm ][ seed ][ argv.mode == "heatmap" ? argv.type : "GMPM" ][ length ];
                    }
                }

                total = total == 0 ? "" : argv.log2 ? Math.log2( total ) : total;
                densities.push( total );
                values.push( "" );

                file_data[ "value"   ][ arm ][ length ] = values;
                file_data[ "density" ][ arm ][ length ] = densities;
            }

            let total = 0;
            for( let i = 0, ilen = total_array.length; i < ilen; i++ ){
                total += total_array[i];

                total_array[i] = total_array[i] == 0 ? "" :
                    argv.log2 ? Math.log2( total_array[i] ) : total_array[i];
            }


            total = total == 0 ? "" : argv.log2 ? Math.log2( total ) : total;
            total_array.push( total );


            file_data[ "value" ][ arm ].push( Array( total_array.length ).fill( "" ));
            file_data[ "value" ][ arm ].reverse();

            file_data[ "density" ][ arm ].push( total_array );
            file_data[ "density" ][ arm ].reverse();
        }

        return file_data;
    }
};
