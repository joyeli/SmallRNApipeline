module.exports = {

    getSeedIndex: function( files, argv ){
        let seed_object = { "5p": {}, "3p":{} };

        for( let file_name in files ){
            for( let arm in files[ file_name ] ){
                for( let seed in files[ file_name ][ arm ] ){
                    if( seed_object[ arm ][ seed ] == null ){ seed_object[ arm ][ seed ] = 0; }
                    for( let length in files[ file_name ][ arm ][ seed ][ "GMPM" ] ){
                        seed_object[ arm ][ seed ] += files[ file_name ][ arm ][ seed ][ "GMPM" ][ length ];
                    }
                }
            }
        }

        let seeds = { "5p": [], "3p":[] };
        let seed_object_array = { "5p": [], "3p":[] };
        let seed_index = { "5p": {}, "3p":{} };

        for( let arm in seed_object ){

            for( let seed in seed_object[ arm ] ){
                seed_object_array[ arm ].push({ Seed: seed, Level: seed_object[ arm ][ seed ]});
                seeds[ arm ].push({ Seed: seed, Level: seed_object[ arm ][ seed ]});
            }

            if( seeds[ arm ].length == 0 ) continue;

            let seed_ref = "";

            if( argv.arm5seq != "" || argv.arm3seq != "" ) {

                if(( argv.arm == "5p3p" && argv.arm5seq != "" && argv.arm3seq == "" ) || ( argv.arm == "5p3p" && argv.arm5seq == "" && argv.arm3seq != "" ) )
                    console.log( "You need to give both \"arm5seq\" and \"arm3seq\" for argv.arm == \"5p3p\"" );

                switch( arm ){
                    case "5p" : seed_ref = argv.arm5seq; break;
                    case "3p" : seed_ref = argv.arm3seq; break;
                }

            }
            else
            {
                seeds[ arm ] = seeds[ arm ].sort( function( a, b ){
                    return a.Level > b.Level ? -1 : 1;
                });

                seed_ref = seeds[ arm ][0].Seed;

                let assemble_parameters = {
                    max_match : 6,
                    min_match : 3,
                    count     : 0,
                };

                seeds[ arm ].shift();
                seed_ref = module.exports.seedAssemble( seed_ref, seeds[ arm ], assemble_parameters );
            }

            for( let i = 0; i <= seed_ref.length -7; ++i ){
                seed_index[ arm ][ seed_ref.substring( i, i+7 )] = i;
            }
        }

        return seed_index;
    },

    seedAssemble: function( seed_ref, seeds, assemble_parameters ){
        if( seeds.length == 0 ){ return seed_ref; }

        let temp_ref  = seeds[0].Seed;
        let max_match = assemble_parameters.max_match;
        let min_match = assemble_parameters.min_match;
        let count     = assemble_parameters.count;

        if( max_match - count < min_match ){
            assemble_parameters.count = 0;
            seeds.shift();

            return module.exports.seedAssemble( seed_ref, seeds, assemble_parameters );
        }

        let temp_ref5 = temp_ref.substring( 0, ( max_match - count ));
        let temp_seed = seed_ref.substring( 0, ( max_match - count ));

        let temp_ref5_idex = seed_ref.substring( seed_ref.length - max_match ).search( temp_ref5 );
        let temp_ref3_idex = temp_ref.substring( temp_ref.length - max_match ).search( temp_seed );

        if( temp_ref5_idex == -1 && temp_ref3_idex == -1 ){
            assemble_parameters.count++;
            seed_ref = module.exports.seedAssemble( seed_ref, seeds, assemble_parameters );
        }
        else if( temp_ref5_idex != -1 ){
            seed_ref = seed_ref.substring( 0, seed_ref.length - max_match + temp_ref5_idex ) + temp_ref;
            seeds.shift();
        }
        else if( temp_ref3_idex != -1 ){
            seed_ref = temp_ref.substring( 0, temp_ref.length - max_match + temp_ref3_idex ) + seed_ref;
            seeds.shift();
        }

        return module.exports.seedAssemble( seed_ref, seeds, assemble_parameters );
    },

    getSeedLabel: function( seed_index, argv ){
        let labelsSeed = { "5p": Array(), "3p": Array() };

        for( let arm in seed_index ){
            if( argv.arm != "5p3p" ){ if( arm != argv.arm ){ continue; }}

            for( let seed in seed_index[ arm ] ){
                labelsSeed[ arm ][ seed_index[ arm ][ seed ]] = seed;
            }

            labelsSeed[ arm ].push( "Total" );
        }

        return labelsSeed;
    }
};
