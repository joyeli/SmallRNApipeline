module.exports = {
    
    log: function( files, datas, argv ){

        let lens = [];
        
        for( let i = argv.minlen; i <= argv.maxlen; ++i ){
            lens.push( i );
        }
        
        let types = Object.keys(
            files[ Object.keys( files )[ 0 ]]
                         [ Object.keys( files  [ Object.keys( files )[ 0 ]] )[ 0 ]]
                         [ Object.keys( files  [ Object.keys( files )[ 0 ]]
                                                       [ Object.keys( files  [ Object.keys( files )[ 0 ]] )[ 0 ]] )[ 0 ]] );
        
        console.log(
                "<table style='border:1px solid black;'><tr><th colspan='2' rowspan='3'"
                + " style='font-family:Sans-serif;'>"
                + ( argv.rename == null ? argv.input : argv.rename ) + "</th>" );
        
        for( let file_name in datas ){
        
            let colspan = 0
        
            for( let arm in datas[ file_name ][ "value" ] ){
                colspan += Object.keys( files[ file_name ][ arm ] ).length;
            }
        
            console.log(
                    "<th colspan='" + colspan
                    + "' style='border-left:1px #666 solid;border-bottom:1px #666 solid;"
                    + "font-family:Sans-serif;font-size:12px;'>"
                    + file_name + "</th>" );
        }
        
        console.log( "</tr><tr>" );
        
        for( let file_name in datas ){
            for( let arm in datas[ file_name ][ "value" ] ){
                if( Object.keys( files[ file_name ][ arm ] ).length == 0 ) continue;
                console.log( "<th colspan='" + Object.keys( files[ file_name ][ arm ] ).length
                        + "' style='font-family:Sans-serif;font-size:12px;border-left:1px #666 solid;'>" + arm + "</th>" );
            }
        }
        
        console.log( "</tr><tr>" );
        
        for( let file_name in datas ){
            for( let arm in datas[ file_name ][ "value" ] ){
        
                let count_arm = 0;
                for( let seed in files[ file_name ][ arm ] ){
        
                    count_arm++;
                    console.log( "<th" + ( count_arm != 1 ? "" : " style='border-left:1px #666 solid;'" )
                            + "><svg width=45 height=85 version='1.1' xmlns='http://www.w3.org/2000/svg'>"
                            + "<text x=-75 y=40 transform='rotate(-70)' font-family='Sans-serif'>"
                            + seed + "</text></svg></th>" );
                }
            }
        }
        
        for( let i = Object.keys( lens ).length -1; i >= 0; --i ){
        
            let count_type = 0;
            console.log(
                    "</tr><tr><th rowspan='8' style='border-top:1px #666 solid;"
                    + "border-right:1px #666 solid;font-family:Sans-serif;"
                    + "'>" + lens[ i ] + "</th>" );
        
            for( let type in types ){
        
                count_type++;
                console.log(
                        "<th style='text-align:right;font-family:Sans-serif;font-size:12px;"
                        + ( count_type != 1 ? "" : "border-top:1px #666 solid;" )
                        + ( types[ type ] != "A_Tail" ? "" : "background-color:#fcc;" )
                        + ( types[ type ] != "C_Tail" ? "" : "background-color:#ccf;" )
                        + ( types[ type ] != "G_Tail" ? "" : "background-color:#ffc;" )
                        + ( types[ type ] != "T_Tail" ? "" : "background-color:#cfc;" )
                        + ( types[ type ] != "Other_Tail" ? "" : "background-color:#eee;" )
                        + "'>" + types[ type ] + "</th>" );
        
                for( let file_name in datas ){
                    for( let arm in datas[ file_name ][ "value" ] ){
        
                        let count_arm = 0;
                        for( let seed in files[ file_name ][ arm ] ){
        
                            count_arm++;
                            console.log(
                                    "<td style='text-align:center;font-family:Sans-serif;font-size:12px;"
                                    + ( count_type != 1 ? "" : "border-top:1px #666 solid;" )
                                    + ( count_arm  != 1 ? "" : "border-left:1px #666 solid;" )
                                    + ( types[ type ] != "A_Tail" ? "" : "background-color:#fcc;" )
                                    + ( types[ type ] != "C_Tail" ? "" : "background-color:#ccf;" )
                                    + ( types[ type ] != "G_Tail" ? "" : "background-color:#ffc;" )
                                    + ( types[ type ] != "T_Tail" ? "" : "background-color:#cfc;" )
                                    + ( types[ type ] != "Other_Tail" ? "" : "background-color:#eee;" )
                                    + "'>"
                                    + ( files[ file_name ][ arm ][ seed ][ types[ type ]][ i ] == 0
                                        ? ""
                                        : ( argv.log2
                                            ? Math.log2( files[ file_name ][ arm ][ seed ][ types[ type ]][ i ] ).toPrecision(3)
                                            : files[ file_name ][ arm ][ seed ][ types[ type ]][ i ] ))
                                    + "</td>"
                                    );
                        }
                    }
                }
        
                console.log( "</tr>" );
            }
        }
        
        console.log( "</tr>" );
        console.log( "</table>" );
    }
};
