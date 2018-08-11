<!DOCTYPE html>
<html>
    <meta charset='utf-8'>
    <script src='https://d3js.org/d3.v5.min.js'></script>
    <script src='https://code.jquery.com/jquery-3.3.1.min.js' ></script>
    <script src='https://code.jquery.com/ui/1.12.1/jquery-ui.js'></script>
    <link rel='stylesheet' href='https://code.jquery.com/ui/1.12.1/themes/smoothness/jquery-ui.css'>
    <link href='https://fonts.googleapis.com/css?family=Source+Code+Pro' rel='stylesheet'>
    <style>.ui-tooltip-content{font-size:15px;font-family:Calibri;}</style>
    <body>
    <?
        Shell_Exec( 'rm /tmp/*' );

        $TSV_File = $_GET['TSV_File'];
        $Annotation_Select = $_GET['Annotation_Select'];

        if( $_GET['TSV_File'] == '' ) $TSV_File = $_POST['TSV_File'];
        if( $_GET['Annotation_Select'] == '' ) $Annotation_Select = $_POST['Annotation_Select'];

        $Length = $_POST['Length'];
        $Filter = $_POST['Filter'];
        $Color_Hight = $_POST['Color_Hight'];
        $Color_Low = $_POST['Color_Low'];
        $Segment_Select = $_POST['Segment_Select'];

        if( $Color_Hight == '' ) $Color_Hight = 'Black';
        if( $Color_Low == '' ) $Color_Low = 'WhiteSmoke';

#<!--================== TSV File ====================-->

        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';

        $TSV = Shell_Exec( 'ls | grep .tsv' );
        $TSV_List = Explode( "\n", $TSV );
        $List_Size = Count( $TSV_List );
        $isRNAfold = false;

        echo '<select name=TSV_File onchange=this.form.submit();>';
        echo '<option '; if($TSV_File=='') echo 'selected'; echo '>Select TSV</option>';

        For( $i = 0; $i < $List_Size-1; ++$i )
        {
            if( $TSV_List[$i] == 'rnafold.tsv' )
            {
                $isRNAfold = true;
                continue;
            }

            echo '<option value='.$TSV_List[$i].' ';

            if( $TSV_File == $TSV_List[$i] ) 
                echo 'selected ';

            echo '>'.$TSV_List[$i].'</option>';
        }
        echo "</select>
            <input type='hidden' name='Length' value='$Length' />
            <input type='hidden' name='Filter' value='$Filter' />
            <input type='hidden' name='Color_Hight' value='$Color_Hight' />
            <input type='hidden' name='Color_Low' value='$Color_Low' />
            <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />
            <input type='hidden' name='Segment_Select' value='$Segment_Select' />
            </form>";

#<!--================== Annotation Select ====================-->

        $inFile = File_get_contents( substr( $TSV_File, 0, strlen( $TSV_File ) -4 ).'.idx' );
        $inFile_Lines = Explode( "\n", $inFile );

        $Total_PPM = 0;
        $PPM_Array = Array();
        $Def_Array = Array();

        $Anno_Array = Array();
        $Data_Array = Array();

        $Ref_Chr    = '';
        $Ref_Start  = '';
        $Ref_Strand = '';
        $Full_Sequc = '';
        $Data_Check = false;

        For( $j = 0; $j < Count( $inFile_Lines )-1; ++$j )
        {
            $anno = Explode( "\t", $inFile_Lines[$j] );
            $Anno_Array[ $anno[0] ] = $anno[1];

            $Anno_Split = Explode( '|', $anno[0] );
            if( Count( $Anno_Split ) != 1 )
            {
                if( Array_Key_Exists( $Anno_Split[0], $Def_Array ))
                {
                    if( $anno[2] > $PPM_Array[ $Def_Array[ $Anno_Split[0] ]] )
                        $Def_Array[ $Anno_Split[0] ] = $anno[0];
                }
                else $Def_Array[ $Anno_Split[0] ] = $anno[0];
            }

            $PPM_Array [ $anno[0] ] = $anno[2];
            if( $Total_PPM < $anno[2] ) $Total_PPM = $anno[2];
        }

        Ksort( $Anno_Array );
        Ksort( $PPM_Array );

        if( Array_Key_Exists( $Annotation_Select, $Def_Array ))
            $Annotation_Select = $Def_Array[ $Annotation_Select ];

        Foreach( $Anno_Array as $anno => $idx )
        {
            if( $Annotation_Select == $anno && $Annotation_Select != '' )
            {
                $TSV = new SplFileObject( $TSV_File );
                $TSV->seek( $idx );

                $annoSeq = Explode( "\t", trim( $TSV->current(), "\n" ));
                $Full_Sequc = $annoSeq[1];
                $Ref_Chr    = $annoSeq[2];
                $Ref_Start  = $annoSeq[3];
                $Ref_Strand = $annoSeq[4];

                while( true )
                {
                    $TSV->next();
                    $tsvData = Explode( "\t", trim( $TSV->current(), "\n" ));
                    if( $tsvData[0] != '' ) break;
                    Array_Push( $Data_Array, $tsvData );
                }

                $TSV = null;
            }
        }

        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';
        echo '<select name=Annotation_Select onchange=this.form.submit();>';
        echo "<option value='' "; if($Annotation_Select=='') echo 'selected'; echo '>Select Annotations</option>';

        Foreach( $Anno_Array as $anno => $idx )
        {
            echo '<option value='.$anno.' ';
            if( $Annotation_Select == $anno )  echo 'selected ';
            echo '>'.$anno.' ('.Number_Format( (float)$PPM_Array[$anno], 2, '.', '' ).'ppm)</option>';
        }

        echo "</select>
            <input type='hidden' name='Length' value='$Length' />
            <input type='hidden' name='Filter' value='$Filter' />
            <input type='hidden' name='Color_Hight' value='$Color_Hight' />
            <input type='hidden' name='Color_Low' value='$Color_Low' />
            <input type='hidden' name='TSV_File' value='$TSV_File' />
            </form>";

        $Total_PPM = 100 / $Total_PPM;
        echo "<script>var select_color_map = d3.scaleLinear().domain([ 0, 100 ]).range([ '$Color_Low', '$Color_Hight' ]);";

        Foreach( $Anno_Array as $anno => $idx )
        {
           echo "$( 'option[value=\"$anno\"]' ).css(
           {
               'background-color': select_color_map( '".($PPM_Array[$anno] * $Total_PPM )."' ),
               'color': '".(( $PPM_Array[$anno] * $Total_PPM ) > 25 ? 'white' : 'black'  )."'
           });";
        }

        echo '</script>';

#<!--=================== RNA-Fold Datas ======================-->

        $RNAfold = Array();

        if( $isRNAfold && $Annotation_Select != '' )
        {
            $inFile = new SplFileObject( 'rnafold.tsv' );

            while( !$inFile->eof() )
            {
                $inFile_Lines = $inFile->fgets();
                $inFile_Line = Explode( "\t", Rtrim( $inFile_Lines ));
                if( $inFile_Line[0] != $Annotation_Select ) continue;

                $Entropy = '';
                $Entropy_Max = 0;

                $Fold = $inFile_Line[5];
                $Sequ = $inFile_Line[4];

                $Fold1 = '';
                $Fold2 = '';
                $Fold3 = '';
                $Fold4 = '';
                $Fold5 = '';

                $Ltmp = '';
                $Mtmp = '';
                $Rtmp = '';

                $L = 0;
                $R = Strlen( $Fold )-1;

                $xEntropy_Array = Array();
                $yEntropy_Array = Array();
                $vEntropy_Array = Array_Slice( $inFile_Line, 6 );

                for( $i = 0; $i < Count( $vEntropy_Array ); ++$i )
                    if( $vEntropy_Array[$i] > $Entropy_Max ) $Entropy_Max = $vEntropy_Array[$i];

                for( $i = 0; $i < Count( $vEntropy_Array ); ++$i )
                {
                    $xEntropy_Array[$i] = 9 + ( $i * 18 );
                    $yEntropy_Array[$i] = 100 - ( $vEntropy_Array[$i] * 100 / $Entropy_Max );
                    $Entropy = $Entropy.$xEntropy_Array[$i].','.$yEntropy_Array[$i].' ';
                }

                for( $i = 0; $i < Strlen( $Fold ); ++$i )
                {
                    if( $Fold[$L] == '(' && $Fold[$R] == ')' )
                    {
                        $Spce = '';
                        $Lspc = '';
                        $Rspc = '';

                        $Lcnt = Strlen( $Ltmp );
                        $Rcnt = Strlen( $Rtmp );
                        $Lnth = $Lcnt < $Rcnt ? $Rcnt : $Lcnt;

                        for( $s = 0; $s < $Lnth; ++$s ) $Spce = $Spce.' ';
                        for( $s = 0; $s < $Rcnt - $Lcnt; ++$s ) $Lspc = $Lspc.' ';
                        for( $s = 0; $s < $Lcnt - $Rcnt; ++$s ) $Rspc = $Rspc.' ';

                        $Fold1 = $Fold1.$Lspc.$Ltmp.' ';
                        $Fold2 = $Fold2.$Spce.$Sequ[$L];
                        $Fold3 = $Fold3.$Spce.'|';
                        $Fold4 = $Fold4.$Spce.$Sequ[$R];
                        $Fold5 = $Fold5.$Rspc.$Rtmp.' ';

                        $Ltmp = '';
                        $Rtmp = '';

                        $L++;
                        $R--;
                    }

                    if( $L == $R )
                    {
                        $Mtmp = $Sequ[$L];
                        break;
                    }

                    if( $Fold[$L] != '(' || $Fold[$R] != ')' )
                    {
                        if( $Fold[$L] != '(' )
                        {
                            $Ltmp = $Ltmp.$Sequ[$L];
                            $L++;
                        }

                        if( $Fold[$R] != ')' )
                        {
                            $Rtmp = $Rtmp.$Sequ[$R];
                            $R--;
                        }
                    }

                    if( $L == ( $R + 1 )) break;
                }

                if( Strlen( $Ltmp ) != Strlen( $Rtmp ))
                {
                    $Lstr = &$Ltmp;
                    $Sstr = &$Rtmp;

                    if( Strlen( $Ltmp ) < Strlen( $Rtmp ))
                    {
                        $Lstr = &$Rtmp;
                        $Sstr = &$Ltmp;
                    }

                    $Dist = Strlen( $Lstr ) - Strlen( $Sstr );

                    if( $Dist % 2 == 0 )
                    {
                        $Sstr = $Sstr.Strrev( Substr( $Lstr, ( Strlen( $Lstr ) - ( $Dist / 2 ))));
                        $Lstr = Substr( $Lstr, 0, Strlen( $Sstr ));
                    }
                    else if( $Mtmp != '' )
                    {
                        $Sstr = $Sstr.$Mtmp;
                        $Dist = Strlen( $Lstr ) - Strlen( $Sstr );
                        $Sstr = $Sstr.Strrev( Substr( $Lstr, ( Strlen( $Lstr ) - ( $Dist / 2 ))));
                        $Lstr = Substr( $Lstr, 0, Strlen( $Sstr ));
                        $Mtmp = '';
                    }
                    else
                    {
                        $Dist = $Dist -1;
                        $Sstr = $Sstr.Strrev( Substr( $Lstr, ( Strlen( $Lstr ) - ( $Dist / 2 ))));
                        $Mtmp = Substr( $Lstr, Strlen( $Sstr ));
                        $Lstr = Substr( $Lstr, 0, Strlen( $Sstr ));
                    }
                }

                for( $i = 0; $i < Strlen( $Ltmp ); ++$i )
                {
                    if( $i < Strlen( $Ltmp ) -1 )
                    {
                        $Fold1 = $Fold1.$Ltmp[$i];
                        $Fold2 = $Fold2.' ';
                        $Fold3 = $Fold3.' ';
                        $Fold4 = $Fold4.' ';
                        $Fold5 = $Fold5.$Rtmp[$i];
                    }
                    else
                    {
                        $Fold1 = $Fold1.' ';
                        $Fold2 = $Fold2.$Ltmp[$i];
                        $Fold3 = $Fold3.( $Mtmp == '' ? ' )' : ( ' '.$Mtmp ));
                        $Fold4 = $Fold4.$Rtmp[$i];
                        $Fold5 = $Fold5.' ';
                    }
                }

                $RNAfold[ 'MFE' ][0] = $inFile_Line[3];
                $RNAfold[ 'Fold0' ][0] = $Fold;
                $RNAfold[ 'Fold1' ][0] = $Fold1;
                $RNAfold[ 'Fold2' ][0] = $Fold2;
                $RNAfold[ 'Fold3' ][0] = $Fold3;
                $RNAfold[ 'Fold4' ][0] = $Fold4;
                $RNAfold[ 'Fold5' ][0] = $Fold5;
                $RNAfold[ 'MaxEntropy' ][0] = $Entropy_Max;
                $RNAfold[ 'Entropy' ][0] = $Entropy;
                $RNAfold[ 'vEntropy' ] = $vEntropy_Array;
                $RNAfold[ 'xEntropy' ] = $xEntropy_Array;
                $RNAfold[ 'yEntropy' ] = $yEntropy_Array;
            }
        }

#<!--=================== Select Segment ======================-->

        $Seg_Start = 0;
        $Seg_End   = 0;

        if( $Full_Sequc != '' && $Data_Array[0][5] == '.' )
        {
            $Segm_Uniqs = Array();
            $Segm_Array = Array();
            $Segm_Range = Array();

            $Segm_Range[0] = 0; // start
            $Segm_Range[1] = 0; // end
            $Segm_Range[2] = 0; // ppm

            $Total_PPM = 0;

            for( $i = 0; $i < Count( $Data_Array ); ++$i )
            {
                $End = $Data_Array[$i][1] + $Data_Array[$i][2]
                    +( $Data_Array[$i][6] != '.' ? Strlen( $Data_Array[$i][6] ) : 0 );

                if( $Data_Array[$i][1] <= $Segm_Range[1] )
                {
                    if( $End > $Segm_Range[1] ) $Segm_Range[1] = $End;
                }
                else
                {
                    if( $Segm_Range[0] != $Segm_Range[1] )
                        Array_push( $Segm_Uniqs, $Segm_Range );

                    for( $j = Count( $Segm_Array ); $j < $i; ++$j )
                        $Segm_Array[$j] = $Segm_Range;

                    $Segm_Range[0] = $Data_Array[$i][1];
                    $Segm_Range[1] = $End;
                    $Segm_Range[2] = 0;
                }

                $Segm_Range[2] += $Data_Array[$i][3];
                if( $Total_PPM < $Data_Array[$i][3] )
                    $Total_PPM = $Data_Array[$i][3];
            }

            $Total_PPM = 100 / $Total_PPM;

            Array_push( $Segm_Uniqs, $Segm_Range );

            if( $Segment_Select == '' )
            {
                $segMax = [ 0, 0, 0 ];
                for( $i = 0; $i < Count( $Segm_Uniqs ); ++$i )
                    if( $segMax[2] < $Segm_Uniqs[$i][2] ) $segMax = $Segm_Uniqs[$i];
                $Segment_Select = $segMax[0].'-'.$segMax[1];
            }

            for( $j = Count( $Segm_Array ); $j < Count( $Data_Array ); ++$j )
                $Segm_Array[$j] = $Segm_Range;

            echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';
            echo '<select name=Segment_Select onchange=this.form.submit();>';
            echo "<option value='' "; if($Segment_Select=='') echo 'selected'; echo '>Select Segments</option>';

            for( $i = 0; $i < Count( $Segm_Uniqs ); ++$i )
            {
                echo '<option value="'.$Segm_Uniqs[$i][0].'-'.$Segm_Uniqs[$i][1].'" ';

                if( $Segment_Select == $Segm_Uniqs[$i][0].'-'.$Segm_Uniqs[$i][1] )
                {
                    echo 'selected ';
                    $Seg_Start = $Segm_Uniqs[$i][0];
                    $Seg_End   = $Segm_Uniqs[$i][1];
                }

                echo '>'
                    .$Ref_Chr.':'
                    .( $Ref_Strand == '+'
                        ? (( (int)$Ref_Start + (int)$Segm_Uniqs[$i][0] ).'-'.( (int)$Ref_Start + (int)$Segm_Uniqs[$i][1] ))
                        : (( (int)$Ref_Start - (int)$Segm_Uniqs[$i][1] ).'-'.( (int)$Ref_Start - (int)$Segm_Uniqs[$i][0] ))
                    ).' ('.Number_Format( (float)$Segm_Uniqs[$i][2], 2, '.', '' ).'ppm)</option>';
            }

            echo "</select>
                <input type='hidden' name='Length' value='$Length' />
                <input type='hidden' name='Filter' value='$Filter' />
                <input type='hidden' name='Color_Hight' value='$Color_Hight' />
                <input type='hidden' name='Color_Low' value='$Color_Low' />
                <input type='hidden' name='TSV_File' value='$TSV_File' />
                <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />
                </form>";

            echo "<script>var select_color_map_seg = d3.scaleLinear().domain([ 0, 100 ]).range([ '$Color_Low', '$Color_Hight' ]);";

            for( $i = 0; $i < Count( $Segm_Uniqs ); ++$i )
                echo "$( 'option[value=\"".$Segm_Uniqs[$i][0].'-'.$Segm_Uniqs[$i][1].' ('.Number_Format( (float)$Segm_Uniqs[$i][2], 2, '.', '' )."ppm)\"]' ).css(
                {
                    'background-color': select_color_map_seg( '".( $Segm_Uniqs[$i][2] * $Total_PPM )."' ),
                    'color': '".(( $Segm_Uniqs[$i][2] * $Total_PPM ) > 25 ? 'white' : 'black'  )."'
                });";

            echo '</script>';
        }

#<!--=================== Specific Length =====================-->

        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';
        echo '<input type=text name=Length size=3 value=';

        if( $Length == '' )
        {
            echo 'Length';
        }
        else echo $Length;

        echo ' onfocus="{this.value=\'\';}" />';

#<!--==================== PPM Filtering ======================-->

        echo '<input type=text name=Filter size=3 value=';

        if( $Filter == '' )
        {
            echo 'PPM';
        }
        else echo $Filter;

        echo ' onfocus="{this.value=\'\';}" />';

#<!--==================== Heatmap Color ======================-->

        echo '<input type=text name=Color_Low size=8 value='.$Color_Low.' onfocus="{this.value=\'\';}" />';
        echo '<input type=text name=Color_Hight size=8 value='.$Color_Hight.' onfocus="{this.value=\'\';}" />';

        echo "
            <input type='hidden' name='TSV_File' value='$TSV_File' />
            <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />
            <input type='hidden' name='Segment_Select' value='$Segment_Select' />
            <input type='submit' value='Submit' /> 
            </form>";

        echo 'MFE: '.$RNAfold[ 'MFE' ][0];

#<!--==================== Make TempFile ======================-->

        $TSV_File_Name = Substr( $TSV_File, 0, ( Strlen( $TSV_File ) -4 ));

        if( $Full_Sequc != '' && ( $Seg_End != 0 || $Data_Array[0][5] != '.' ))
        {
            $Temp = Tempnam( '/tmp', $TSV_File.'_'.$Annotation_Select.'_len'.$Length.'_ppm'.$Filter.'_JSON_' );
            $Ftemp = Fopen( $Temp, 'w' );

            Fwrite( $Ftemp, "{\n" );
            Fwrite( $Ftemp, "  \"Annotation\" : \"".$Annotation_Select."\",\n" );
            Fwrite( $Ftemp, "  \"Sequence\"   : \"".
                ( $Seg_End == 0 ? $Full_Sequc :
                    ( $Data_Array[0][5] == '.'
                        ? Substr( $Full_Sequc, $Seg_Start - 10, ( $Seg_End - $Seg_Start + 20 ))
                        : Substr( $Full_Sequc, $Seg_Start, ( $Seg_End - $Seg_Start ))
                    )
                )."\",\n" );

            if( $isRNAfold )
            {
                Fwrite( $Ftemp, "  \"Fold0\"      : \"".$RNAfold[ 'Fold0' ][0]."\",\n" );
                Fwrite( $Ftemp, "  \"Fold1\"      : \"".$RNAfold[ 'Fold1' ][0]."\",\n" );
                Fwrite( $Ftemp, "  \"Fold2\"      : \"".$RNAfold[ 'Fold2' ][0]."\",\n" );
                Fwrite( $Ftemp, "  \"Fold3\"      : \"".$RNAfold[ 'Fold3' ][0]."\",\n" );
                Fwrite( $Ftemp, "  \"Fold4\"      : \"".$RNAfold[ 'Fold4' ][0]."\",\n" );
                Fwrite( $Ftemp, "  \"Fold5\"      : \"".$RNAfold[ 'Fold5' ][0]."\",\n" );
                Fwrite( $Ftemp, "  \"MFE\"        : ".$RNAfold[ 'MFE' ][0].",\n" );
                Fwrite( $Ftemp, "  \"MaxEntropy\" : ".$RNAfold[ 'MaxEntropy' ][0].",\n" );
                Fwrite( $Ftemp, "  \"Entropy\"    : \"".$RNAfold[ 'Entropy' ][0]."\",\n" );

                Fwrite( $Ftemp, "  \"vEntropy\"   : [" );
                for( $i = 0; $i < Count( $RNAfold[ 'vEntropy' ] ); ++$i )
                    Fwrite( $Ftemp, $RNAfold[ 'vEntropy' ][$i].( $i < Count( $RNAfold[ 'vEntropy' ] )-1 ? ',' : "],\n" ));

                Fwrite( $Ftemp, "  \"xEntropy\"   : [" );
                for( $i = 0; $i < Count( $RNAfold[ 'xEntropy' ] ); ++$i )
                    Fwrite( $Ftemp, $RNAfold[ 'xEntropy' ][$i].( $i < Count( $RNAfold[ 'xEntropy' ] )-1 ? ',' : "],\n" ));

                Fwrite( $Ftemp, "  \"yEntropy\"   : [" );
                for( $i = 0; $i < Count( $RNAfold[ 'yEntropy' ] ); ++$i )
                    Fwrite( $Ftemp, $RNAfold[ 'yEntropy' ][$i].( $i < Count( $RNAfold[ 'yEntropy' ] )-1 ? ',' : "],\n" ));
            }

            Fwrite( $Ftemp, "  \"RefChr\"     : \"".$Ref_Chr."\",\n" );
            Fwrite( $Ftemp, "  \"RefStart\"   : ".$Ref_Start.",\n" );
            Fwrite( $Ftemp, "  \"RefStrand\"  : \"".$Ref_Strand."\",\n" );

            Fwrite( $Ftemp, "  \"Reads\"      : [\n" );

            $isFirst = true;

            for( $i = 0; $i < Count( $Data_Array ); ++$i )
            {
                if( $Seg_End != 0 && $Segm_Array[$i][0] != $Seg_Start && $Segm_Array[$i][1] != $Seg_End )
                    continue;

                if( !$isFirst ) Fwrite( $Ftemp, ",\n" );

                Fwrite( $Ftemp, "    {\n" );
                Fwrite( $Ftemp, "      \"Index\"    : "  .( $Data_Array[$i][1] - $Seg_Start + ( $Data_Array[0][5] == '.' ? 10 : 0 )).  ",\n" );
                Fwrite( $Ftemp, "      \"Length\"   : "  .$Data_Array[$i][2].  ",\n" );
                Fwrite( $Ftemp, "      \"PPM\"      : "  .$Data_Array[$i][3].  ",\n" );
                Fwrite( $Ftemp, "      \"isRMSK\"   : \"".$Data_Array[$i][4]."\",\n" );
                Fwrite( $Ftemp, "      \"Arm\"      : \"".$Data_Array[$i][5]."\",\n" );
                Fwrite( $Ftemp, "      \"SegStart\" : \"".  $Segm_Array[$i][0]  . "\",\n" );
                Fwrite( $Ftemp, "      \"SegEnd\"   : \"".  $Segm_Array[$i][1]  . "\",\n" );
                Fwrite( $Ftemp, "      \"Tail\"     : \"".( $Data_Array[$i][6] != "." ? $Data_Array[$i][6] : "" )."\",\n" );
                Fwrite( $Ftemp, "      \"MDtag\"    : \"".( $Data_Array[$i][7] != "." ? $Data_Array[$i][7] : "" )."\",\n" );
                Fwrite( $Ftemp, "      \"TCtag\"    : \"".( $Data_Array[$i][8] != "." ? $Data_Array[$i][8] : "" )."\",\n" );

                if( ( $Length != '' && $Length != 'Length' && $Length != $Data_Array[$i][2] ) ||
                    ( $Filter != '' && $Filter != 'PPM'    && $Filter >  $Data_Array[$i][3] ) )
                    Fwrite( $Ftemp, "      \"Filter\" : \"Y\"\n" );
                else
                    Fwrite( $Ftemp, "      \"Filter\" : \"N\"\n" );

                Fwrite( $Ftemp, "    }" );

                $isFirst = false;
            }

            Fwrite( $Ftemp, "\n  ]\n" );
            Fwrite( $Ftemp, "}\n" );
            Fclose( $Ftemp );

#<!--==================== Sq Align Plot ======================-->

            echo "<script>

                var color_map = d3.scaleLinear().domain([ 0, 100 ]).range([ '$Color_Low', '$Color_Hight' ]);

                var datas = Array();
                var exprs = Array();
                var etpys = Array();
                var seeds = Array();

                var anno = '';
                var sequ = '';

                var colorA = 'red';
                var colorC = 'blue';
                var colorG = 'goldenrod';
                var colorT = 'green';
                var colorP = 'darkgray';

                var fold_height = 200;
                var exp_height = 120;
                var seq_height = 50;
                var read_height = 20;
                var spc_num = 18;
                var expr_max = 0;
                var shift_top1 = 30 + fold_height;
                var shift_top2 = 38;

                var seg_start = $Seg_Start;
                var seg_end   = $Seg_End;

                $.getJSON( '$Temp', function( json )
                {
                    anno = json[ 'Annotation' ];
                    sequ = json[ 'Sequence'   ];

                    for( i in json[ 'vEntropy' ] )
                        etpys[i] = json[ 'vEntropy' ][i];

                    var expr_array = Array();
                    var heat_array = Array();

                    var heat_max = 0;
                    var lb_width = 1;

                    var start = 0;
                    var end   = 0;
                    var extra = 0;

                    $.each( json[ 'Reads' ], function( idx, read )
                    {
                        end = read[ 'Index' ] + read[ 'Length' ] + read[ 'Tail' ].length;

                        if(( end - json[ 'Sequence' ].length ) > 0 &&
                           ( end - json[ 'Sequence' ].length ) > extra )
                        {
                            extra = end - json[ 'Sequence' ].length;
                        }

                        if( !(( read[ 'Index' ] + 4 ) in seeds ))
                        {
                            seeds[ read[ 'Index' ] + 4 ] = Array();
                            seeds[ read[ 'Index' ] + 4 ][ 'Name' ]
                                = anno + ".( $Data_Array[0][5] != '.' ? ( "'-' + read[ 'Arm' ]" ) : "''" )." + '_' + sequ.substr(( read[ 'Index' ] + 1 ), 7 );

                            seeds[ read[ 'Index' ] + 4 ][ 'GMPM' ] = 0; 
                            seeds[ read[ 'Index' ] + 4 ][  'GM'  ] = 0;
                            seeds[ read[ 'Index' ] + 4 ][  'PM'  ] = 0;
                            seeds[ read[ 'Index' ] + 4 ][ 'Counts' ] = 0;
                            seeds[ read[ 'Index' ] + 4 ][ 'Length' ] = 0;
                        }

                        seeds[ read[ 'Index' ] + 4 ][ 'GMPM' ] += read[ 'PPM' ];
                        seeds[ read[ 'Index' ] + 4 ][  'GM'  ] += ( read[ 'Tail' ] == '' ? read[ 'PPM' ] : 0 );
                        seeds[ read[ 'Index' ] + 4 ][  'PM'  ] += ( read[ 'Tail' ] != '' ? read[ 'PPM' ] : 0 );
                        seeds[ read[ 'Index' ] + 4 ][ 'Counts' ]++;
                        seeds[ read[ 'Index' ] + 4 ][ 'Length' ]
                            = (( read[ 'Length' ] + read[ 'Tail' ].length ) > seeds[ read[ 'Index' ] + 4 ][ 'Length' ]
                            ?  ( read[ 'Length' ] + read[ 'Tail' ].length )
                            :    seeds[ read[ 'Index' ] + 4 ][ 'Length' ] )
                            ;
                    });

                    for( var i = 0; i < json[ 'Sequence' ].length + extra; ++i )
                    {
                        expr_array[ i ] = Array();
                        expr_array[ i ][ 'PPM' ] = 0;
                        expr_array[ i ][  'P'  ] = 0;
                        expr_array[ i ][  'A'  ] = 0;
                        expr_array[ i ][  'C'  ] = 0;
                        expr_array[ i ][  'G'  ] = 0;
                        expr_array[ i ][  'T'  ] = 0;
                    }

                    $.each( json[ 'Reads' ], function( idx, read )
                    {
                        if( read[ 'Filter' ] == 'Y' ) return;

                        datas[ idx ] = read;
                        heat_array[ idx ] = read[ 'PPM' ];

                        start = read[ 'Index' ];
                        end   = read[ 'Index' ] + read[ 'Length' ];

                        for( var i = start; i < ( end + read[ 'Tail' ].length ); ++i )
                        {
                            expr_array[ i ][ 'PPM' ] += read[ 'PPM' ];

                            if( i < end ) expr_array[ i ][  'P'  ] += read[ 'PPM' ];
                        }

                        for( var i = 0; i < read[ 'Tail' ].length; ++i )
                        {
                            if( read[ 'Tail' ].charAt(i) == 'A' ) expr_array[ i + end ][ 'A' ] += read[ 'PPM' ];
                            if( read[ 'Tail' ].charAt(i) == 'C' ) expr_array[ i + end ][ 'C' ] += read[ 'PPM' ];
                            if( read[ 'Tail' ].charAt(i) == 'G' ) expr_array[ i + end ][ 'G' ] += read[ 'PPM' ];
                            if( read[ 'Tail' ].charAt(i) == 'T' ) expr_array[ i + end ][ 'T' ] += read[ 'PPM' ];
                        }
                    });

                    $.each( expr_array, function( idx )
                    {
                        if( expr_max < expr_array[ idx ][ 'PPM' ] ) expr_max = expr_array[ idx ][ 'PPM' ];
                    });

                    $.each( expr_array, function( idx )
                    {
                        expr_array[ idx ][ 'PPM' ] = expr_array[ idx ][ 'PPM' ] / expr_max * 100;

                        var total = expr_array[ idx ][ 'P' ]
                                  + expr_array[ idx ][ 'A' ]
                                  + expr_array[ idx ][ 'C' ]
                                  + expr_array[ idx ][ 'G' ]
                                  + expr_array[ idx ][ 'T' ];

                        if( total != 0 )
                        {
                            expr_array[ idx ][ 'P' ] = expr_array[ idx ][ 'P' ] / total;
                            expr_array[ idx ][ 'A' ] = expr_array[ idx ][ 'A' ] / total;
                            expr_array[ idx ][ 'C' ] = expr_array[ idx ][ 'C' ] / total;
                            expr_array[ idx ][ 'G' ] = expr_array[ idx ][ 'G' ] / total;
                            expr_array[ idx ][ 'T' ] = expr_array[ idx ][ 'T' ] / total;
                        }
                    });

                    $.each( heat_array, function( idx )
                    {
                        if( heat_max < heat_array[ idx ] )
                            heat_max = heat_array[ idx ];
                    });

                    $.each( heat_array, function( idx )
                    {
                        heat_array[ idx ] = heat_array[ idx ] / heat_max * 100;
                        if( heat_array[ idx ] == 0 ) heat_array[ idx ] = 0.01;
                    });

                    $( 'body' ).append( \"<div id='chart'></div>\" );
                    $( '#chart' ).css({
                        'margin-top': '10px',
                        });

                    $( '#chart' ).append( \"<div id='foldchart'></div>\" );
                    $( '#foldchart' ).css({
                        'height': fold_height,
                        'width': (( json[ 'Sequence' ].length + extra + lb_width )* spc_num ) + 'px',
                        });

                    var foldidx = 0;
                    var foldidxs = [];

                    for( var f = 1; f <= 5; f++ )
                        foldidxs[f] = [];

                    for( var i = 0; i < json[ 'Fold1' ].length; ++i )
                    {
                        for( var f = 1; f <= 2; f++ )
                        {
                            if( json[ 'Fold' + f ].charAt(i) != ' ' )
                            {
                                foldidx++;
                                foldidxs[f][i] = foldidx;
                            }
                            else foldidxs[f][i] = 0;
                        }
                    }

                    for( var i = 0; i < json[ 'Fold3' ].length; ++i )
                    {
                        if( json[ 'Fold3' ].charAt(i) != ' ' && json[ 'Fold3' ].charAt(i) != '|' && json[ 'Fold3' ].charAt(i) != ')' )
                        {
                            foldidx++;
                            foldidxs[3][i] = foldidx;
                        }
                        else foldidxs[3][i] = 0;
                    }

                    for( var i = json[ 'Fold5' ].length -1; i >= 0; --i )
                    {
                        for( var f = 4; f <= 5; f++ )
                        {
                            if( json[ 'Fold' + f ].charAt(i) != ' ' )
                            {
                                foldidx++;
                                foldidxs[f][i] = foldidx;
                            }
                            else foldidxs[f][i] = 0;
                        }
                    }

                    for( var f = 1; f <= 5; f++ )
                    {
                        var foldx = 'Fold' + f;
                        var left_shit = (( json[ 'Sequence' ].length * spc_num ) /2 ) - ( json[ foldx ].length * 18 / 2 );

                        $( '#foldchart' ).append( \"<div id='\" + foldx + \"'></div>\" );
                        $( '#' + foldx ).css({
                            'position': 'relative',
                            'left': left_shit + 'px',
                            'width': ( json[ 'Sequence' ].length * spc_num ) - left_shit + 'px',
                            'font-family': 'Source Code Pro, monospace',
                            'font-size': '30px',
                            });

                        for( var i = 0; i < json[ foldx ].length; ++i )
                        {
                            $( '#' + foldx ).append( \"<div id='fold\" + ( foldidxs[f][i] == 0 ? ( f*1000+i ) : foldidxs[f][i] ) + \"' title='fold\" + foldidxs[f][i] + \"'>\" + json[ foldx ].charAt(i) + '</div>' );
                            $( '#fold' + ( foldidxs[f][i] == 0 ? ( f*1000+i ) : foldidxs[f][i] )).css({ 'white-space': 'pre' });
                            if( i != json[ foldx ].length-1 )
                                $( '#fold' + ( foldidxs[f][i] == 0 ? ( f*1000+i ) : foldidxs[f][i] )).css({ 'float': 'left' });
                        }
                    }

                    $( '#chart' ).append( \"<div id='expression'></div>\" );
                    $( '#expression' ).css({
                        'height': '120px',
                        'width': (( json[ 'Sequence' ].length + extra + lb_width )* spc_num ) + 'px',
                        });

                    $( '#expression' ).append( \"<div id='exprchart'></div>\" );
                    $( '#exprchart' ).css({
                        'height': '104px',
                        'width': (( json[ 'Sequence' ].length + extra + lb_width )* spc_num ) + 'px',
                        });

                    $( '#expression' ).append( $(\"<svg id='etpchart' xmlns='http://www.w3.org/2000/svg'><polyline id='entropy' points='\" + json[ 'Entropy' ] + \"'/></svg>\"))
                    $( '#etpchart' ).css({
                        'position': 'absolute',
                        'z-index': '-1',
                        'height': '100px',
                        'width': (( json[ 'Sequence' ].length + extra + lb_width )* spc_num ) + 'px',
                        'top': 7 + shift_top1 + 'px',
                        });

                    $( '#entropy' ).css({
                        'fill': 'none',
                        'stroke': 'Paleturquoise',
                        'stroke-width': 2
                        });

                    $.each( expr_array, function( idx, expr )
                    {
                        exprs[ idx ] = expr;

                        $( '#exprchart' ).append( \"<div id='exprtitle\" + idx + \"' title='expr\" + idx + \"'></div>\" );
                        $( '#exprtitle' + idx ).css({
                            'height': '100px',
                            'width': spc_num + 'px',
                            'float': 'left'
                            });

                        $( '#exprtitle' + idx ).append( \"<div id='expr\" + idx + \"'></div>\" );
                        $( '#expr' + idx ).css({
                            'margin-top': ( 100 - expr[ 'PPM' ] ) + 'px',
                            'width': spc_num + 'px',
                            'position': 'relative',
                            'float': 'left',
                            'z-index': '-1'
                            });

                        $( '#expr' + idx ).append( \"<div id='expr\" + idx + \"T'></div>\" );
                        $( '#expr' + idx ).append( \"<div id='expr\" + idx + \"G'></div>\" );
                        $( '#expr' + idx ).append( \"<div id='expr\" + idx + \"C'></div>\" );
                        $( '#expr' + idx ).append( \"<div id='expr\" + idx + \"A'></div>\" );
                        $( '#expr' + idx ).append( \"<div id='expr\" + idx + \"P'></div>\" );

                        $( '#expr' + idx + 'T' ).css({ 'height': ( expr[ 'PPM' ] * expr[ 'T' ]) + 'px', 'background': colorT });
                        $( '#expr' + idx + 'G' ).css({ 'height': ( expr[ 'PPM' ] * expr[ 'G' ]) + 'px', 'background': colorG });
                        $( '#expr' + idx + 'C' ).css({ 'height': ( expr[ 'PPM' ] * expr[ 'C' ]) + 'px', 'background': colorC });
                        $( '#expr' + idx + 'A' ).css({ 'height': ( expr[ 'PPM' ] * expr[ 'A' ]) + 'px', 'background': colorA });
                        $( '#expr' + idx + 'P' ).css({ 'height': ( expr[ 'PPM' ] * expr[ 'P' ]) + 'px', 'background': colorP });
                    });

                    $( '#exprchart' ).append( \"<div id='exprlabel'></div>\" );
                    $( '#exprlabel' ).css({
                        'width': '1px',
                        'height': '100px',
                        'background': 'black',
                        'float': 'left'
                        });

                    $( '#exprlabel' ).append( \"<div id='exprlabeltop'>-\" + json[ 'MaxEntropy' ].toFixed(2) + '</div>' );
                    $( '#exprlabeltop' ).css({
                        'font-family': 'Source Code Pro, monospace',
                        'font-size': '15px',
                        'position': 'absolute',
                        'top': 0 + shift_top1 + 'px',
                        'float': 'left'
                        });

                    $( '#exprlabel' ).append( \"<div id='exprlabeldown'>-0.00</div>\" );
                    $( '#exprlabeldown' ).css({
                        'font-family': 'Source Code Pro, monospace',
                        'font-size': '15px',
                        'position': 'absolute',
                        'top': 100 + shift_top1 + 'px',
                        'float': 'left'
                        });

                    $( '#expression' ).append( \"<div id='basecount1'></div>\" );
                    $( '#basecount1' ).css({ 'height': '10px' });

                    for( var i = 1; i < json[ 'Sequence' ].length + extra; ++i ) if( i % 5 == 0 )
                    {
                        $( '#basecount1' ).append( \"<div id='base1-\" + i + \"' align='right'>\" + ( i == 5 ? '05' : i ) + '</div>' );
                        $( '#base1-' + i ).css({
                            'width': ( 5 * spc_num ) + 'px',
                            'font-family': 'Source Code Pro, monospace',
                            'font-size': '10px',
                            'float': 'left'
                            });
                    }

                    $( '#chart' ).append( \"<div id='sequence'></div>\" );
                    $( '#sequence' ).css({
                        'width': (( json[ 'Sequence' ].length + extra + lb_width )* spc_num ) + 'px',
                        'height': seq_height + 'px',
                        });

                    $( '#sequence' ).append( \"<div id='selectbox'></div>\" );
                    $( '#selectbox' ).css({
                        'border': '5px solid gold',
                        'border-radius': '12px',
                        'position': 'absolute',
                        'pointer-events': 'none',
                        'height': '68px',
                        'display': 'none'
                        });

                    $( '#sequence' ).append( \"<div id='seqs'></div>\" );

                    $( '#seqs' ).css({
                        'font-family': 'Source Code Pro, monospace',
                        'font-size': '30px',
                        });

                    for( var i = 0; i < json[ 'Sequence' ].length; ++i )
                    {
                        if( i in seeds )
                        {
                            var MDarray = Array();
                            var MDseed = '';
                            var RMSK = '';

                            $.each( json[ 'Reads' ], function( idx, read )
                            {
                                if( read[ 'Index' ] != i ) return;
                                RMSK = read[ 'isRMSK' ];
                                var pos = '';

                                for( var j = 0; j < read[ 'MDtag' ].length; ++j )
                                {
                                    if( !$.isNumeric( read[ 'MDtag' ][j] ))
                                    {
                                        MDarray[ pos ] = read[ 'MDtag' ][j];
                                        pos = '';
                                    }
                                    else pos += read[ 'MDtag' ][j];
                                }
                            });

                            for( var j = 1; j < 8; ++j )
                                if( j in MDarray ) MDseed += ( j -1 ) + MDarray[j];

                            var anno_split = seeds[i][ 'Name' ].split( '|' );
                            var annotation = anno_split[0];

                            if( anno_split.length != 1 )
                            {
                                anno_split = anno_split[1].split( '_' );
                                annotation = annotation + '_' + anno_split[1]
                            }

                            annotation = annotation
                                + ( MDseed == '' ? '' : ( '|' + MDseed ))
                                + ( RMSK == 'N' ? '' : '!' );

                            $( '#seqs' ).append( \"<a id='seq\" + i + \"linked' target='_blank' href='../LenDist/index.php?TSV_File=$TSV_File_Name&annotation_select=\" + annotation + \"'></a>\" );
                            $( '#seq' + i + 'linked' ).css({
                                'text-decoration': 'none',
                                'color': 'black'
                                });

                            $( '#seq' + i + 'linked' ).append( \"<div id='seq\" + i + \"' title='seqs\" + ( i + 1 ) + \"'>\" + json[ 'Sequence' ].charAt(i) + '</div>' );
                        }
                        else $( '#seqs' ).append( \"<div id='seq\" + i + \"' title='seqs\" + ( i + 1 ) + \"'>\" + json[ 'Sequence' ].charAt(i) + '</div>' );

                        if( i != json[ 'Sequence' ].length-1 )
                            $( '#seq' + i ).css({ 'float': 'left' });
                    }

                    $( '#sequence' ).append( \"<div id='folds'></div>\" );

                    $( '#folds' ).css({
                        'font-family': 'Source Code Pro, monospace',
                        'font-size': '30px',
                        });

                    for( var i = 0; i < json[ 'Fold0' ].length; ++i )
                    {
                        $( '#folds' ).append( \"<div id='fold0\" + i + \"' title='seqs\" + ( i + 1 ) + \"'>\" + json[ 'Fold0' ].charAt(i) + '</div>' );
                        if( i != json[ 'Fold0' ].length-1 )
                            $( '#fold0' + i ).css({ 'float': 'left' });
                    }

                    $( '#sequence' ).append( \"<div id='basecount2'></div>\" );
                    $( '#basecount2' ).css({ 'height': '10px' });

                    for( var i = 1; i < json[ 'Sequence' ].length + extra; ++i ) if( i % 5 == 0 )
                    {
                        $( '#basecount2' ).append( \"<div id='base3-\" + i + \"' align='right'>\" + ( i == 5 ? '05' : i ) + '</div>' );
                        $( '#base3-' + i ).css({
                            'width': ( 5 * spc_num ) + 'px',
                            'font-family': 'Source Code Pro, monospace',
                            'font-size': '10px',
                            'float': 'left'
                            });
                    }

                    $( '#chart' ).append( \"<div id='seqalign'></div>\" );
                    $( '#seqalign' ).css({
                        'width': (( json[ 'Sequence' ].length + extra + lb_width )* spc_num ) + 'px',
                        'float': 'left'
                        });

                    $( '#seqalign' ).append( \"<div id='reads'></div>\" );
                    $( '#reads' ).css({
                        'width': (( json[ 'Sequence' ].length + extra )* spc_num ) + 'px',
                        'float': 'left'
                        });

                    $( '#seqalign' ).append( \"<div id='readslabel'></div>\" );
                    $( '#readslabel' ).css({
                        'width': '5px',
                        'height': ( datas.filter( Boolean ).length * read_height ) + 'px',
                        'background': 'linear-gradient( to top, $Color_Low 0%, $Color_Hight 100% )',
                        'float': 'left'
                        });

                    $( '#readslabel' ).append( \"<div id='readslabeltop'>-\" + heat_max.toFixed(2) + '</div>' );
                    $( '#readslabeltop' ).css({
                        'font-family': 'Source Code Pro, monospace',
                        'font-size': '15px',
                        'position': 'absolute',
                        'top': seq_height + exp_height + shift_top1 + shift_top2 + 'px',
                        'float': 'left'
                        });

                    if( datas.filter( Boolean ).length > 4 )
                    {
                        $( '#readslabel' ).append( \"<div id='readslabelq2'>-\" + ( heat_max * 0.75 ).toFixed(2) + '</div>' );
                        $( '#readslabelq2' ).css({
                            'font-family': 'Source Code Pro, monospace',
                            'font-size': '15px',
                            'position': 'absolute',
                            'top': ( datas.filter( Boolean ).length * read_height * 0.25 ) + seq_height + exp_height + shift_top1 + shift_top2 + 'px',
                            'float': 'left'
                            });
                    }

                    if( datas.filter( Boolean ).length > 1 )
                    {
                        $( '#readslabel' ).append( \"<div id='readslabelmid'>-\" + ( heat_max * 0.5 ).toFixed(2) + '</div>' );
                        $( '#readslabelmid' ).css({
                            'font-family': 'Source Code Pro, monospace',
                            'font-size': '15px',
                            'position': 'absolute',
                            'top': ( datas.filter( Boolean ).length * read_height * 0.5 ) + seq_height + exp_height + shift_top1 + shift_top2 + 'px',
                            'float': 'left'
                            });
                    }

                    if( datas.filter( Boolean ).length > 4 )
                    {
                        $( '#readslabel' ).append( \"<div id='readslabelq1'>-\" + ( heat_max * 0.25 ).toFixed(2) + '</div>' );
                        $( '#readslabelq1' ).css({
                            'font-family': 'Source Code Pro, monospace',
                            'font-size': '15px',
                            'position': 'absolute',
                            'top': ( datas.filter( Boolean ).length * read_height * 0.75 ) + seq_height + exp_height + shift_top1 + shift_top2 + 'px',
                            'float': 'left'
                            });
                    }

                    $( '#readslabel' ).append( \"<div id='readslabeldown'>-0.00</div>\" );
                    $( '#readslabeldown' ).css({
                        'font-family': 'Source Code Pro, monospace',
                        'font-size': '15px',
                        'position': 'absolute',
                        'top': ( datas.filter( Boolean ).length * read_height ) + seq_height + exp_height + shift_top1 + shift_top2 + 'px',
                        'float': 'left'
                        });

                    var last_idx = -1;
                    var reads_count = 0;

                    $.each( json[ 'Reads' ], function( idx, read )
                    {
                        if( read[ 'Filter' ] == 'Y' ) return;
                        reads_count++;

                        var pos = '';
                        var MDseed = '';
                        var MDarray = Array();
                        var TCarray = Array();

                        for( var i = 0; i < read[ 'MDtag' ].length; ++i )
                        {
                            if( !$.isNumeric( read[ 'MDtag' ][i] ))
                            {
                                MDarray[ pos ] = read[ 'MDtag' ][i];
                                pos = '';
                            }
                            else pos += read[ 'MDtag' ][i];
                        }

                        for( var i = 0; i < read[ 'TCtag' ].length; ++i )
                        {
                            if( !$.isNumeric( read[ 'TCtag' ][i] ))
                            {
                                TCarray[ pos ] = read[ 'TCtag' ][i];
                                pos = '';
                            }
                            else pos += read[ 'TCtag' ][i];
                        }

                        for( var i = 1; i < 8; ++i )
                            if( i in MDarray ) MDseed += ( i -1 ) + MDarray[i];

                        if( last_idx != read[ 'Index' ] )
                        {
                            last_idx = read[ 'Index' ];
                            $( '#reads' ).append( \"<div id='selectedseed\" + last_idx + \"'></div>\" );

                            $( '#selectedseed' + last_idx ).css({
                                'position': 'relative',
                                'left': last_idx * spc_num + 'px',
                                'width': seeds[ last_idx + 4 ][ 'Length' ] * spc_num + 'px',
                                });

                            var anno_split = seeds[ last_idx + 4 ][ 'Name' ].split( '|' );
                            var annotation = anno_split[0];

                            if( anno_split.length != 1 )
                            {
                                anno_split = anno_split[1].split( '_' );
                                annotation = annotation + '_' + anno_split[1]
                            }

                            $( '#selectedseed' + last_idx ).append( \"<a id='read\" + last_idx + \"linked' target='_blank' href='../LenDist/index.php?TSV_File=$TSV_File_Name&annotation_select=\"
                                + annotation
                                + ( MDseed == '' ? '' : ( '|' + MDseed ))
                                + ( read[ 'isRMSK' ] == 'N' ? '' : '!' )
                                + \"'></a>\" );

                            $( '#read' + last_idx + 'linked' ).css({ 'text-decoration': 'none' });
                        }

                        $( '#read' + last_idx + 'linked' ).append( \"<div id='read\" + idx + \"' title=Read\" + idx + \"></div>\" );
                        $( '#read' + idx ).css({
                            'width': (( read[ 'Length' ] + read[ 'Tail' ].length ) * spc_num ) + 'px',
                            });

                        pos = 0;

                        for( var i = 0; i < read[ 'Length' ]; ++i )
                        {
                            if( i in MDarray )
                            {
                                $( '#read' + idx ).append( \"<div id='read\" + idx + 'gm' + i + \"'></div>\" );
                                $( '#read' + idx + 'gm' + i ).css({
                                    'width': ( pos * spc_num ) + 'px',
                                    'height': read_height + 'px',
                                    'background': color_map( heat_array[ idx ]),
                                    'position': 'relative',
                                    'float': 'left',
                                    'z-index': '-1'
                                    });

                                pos = 0;

                                $( '#read' + idx ).append( \"<div id='read\" + idx + 'md' + i + \"'>\" + MDarray[i] + '</div>' );
                                $( '#read' + idx + 'md' + i ).css({
                                    'font-family': 'Source Code Pro, monospace',
                                    'font-weight': '100',
                                    'font-size': '18px',
                                    'text-align': 'center',
                                    'height': read_height + 'px',
                                    'width': spc_num + 'px',
                                    'background': color_map( heat_array[ idx ]),
                                    'color': 'lightpink',
                                    'position': 'relative',
                                    'float': 'left',
                                    'z-index': '-1'
                                    });
                            }
                            else if( i in TCarray )
                            {
                                $( '#read' + idx ).append( \"<div id='read\" + idx + 'gm' + i + \"'></div>\" );
                                $( '#read' + idx + 'gm' + i ).css({
                                    'width': ( pos * spc_num ) + 'px',
                                    'height': read_height + 'px',
                                    'background': color_map( heat_array[ idx ]),
                                    'position': 'relative',
                                    'float': 'left',
                                    'z-index': '-1'
                                    });

                                pos = 0;

                                $( '#read' + idx ).append( \"<div id='read\" + idx + 'tc' + i + \"'>\" + TCarray[i] + '</div>' );
                                $( '#read' + idx + 'tc' + i ).css({
                                    'font-family': 'Source Code Pro, monospace',
                                    'font-weight': '100',
                                    'font-size': '18px',
                                    'text-align': 'center',
                                    'height': read_height + 'px',
                                    'width': spc_num + 'px',
                                    'background': color_map( heat_array[ idx ]),
                                    'color': 'lightblue',
                                    'position': 'relative',
                                    'float': 'left',
                                    'z-index': '-1'
                                    });
                            }
                            else
                            {
                                pos++;

                                if( i == read[ 'Length' ] -1 )
                                {
                                    $( '#read' + idx ).append( \"<div id='read\" + idx + 'gm' + i + \"'></div>\" );
                                    $( '#read' + idx + 'gm' + i ).css({
                                        'width': ( pos * spc_num ) + 'px',
                                        'height': read_height + 'px',
                                        'background': color_map( heat_array[ idx ]),
                                        'position': 'relative',
                                        'float': 'left',
                                        'z-index': '-1'
                                        });

                                    if( read[ 'Tail' ] == '' )
                                    {
                                        $( '#read' + idx ).append( \"<div id='read\" + idx + 'gmEnd' + \"'></div>\" );
                                        $( '#read' + idx + 'gmEnd' ).css({
                                            'width': ( 0 * spc_num ) + 'px',
                                            'height': read_height + 'px',
                                            'background': color_map( heat_array[ idx ]),
                                            'position': 'relative',
                                            'z-index': '-1'
                                            });
                                    }
                                }
                            }
                        }

                        for( var i = 0; i < read[ 'Tail' ].length; ++i )
                        {
                            $( '#read' + idx ).append( \"<div id='read\" + idx + 'pm' + i + \"'>\" + read[ 'Tail' ].charAt(i) + '</div>' );
                            $( '#read' + idx + 'pm' + i ).css({
                                'font-family': 'Source Code Pro, monospace',
                                'font-size': '18px',
                                'text-align': 'center',
                                'height': read_height + 'px',
                                'position': 'relative',
                                'z-index': '-1'
                                });

                            if( read[ 'Tail' ].charAt(i) == 'A' ) $( '#read' + idx + 'pm' + i ).css({ 'color': colorA });
                            if( read[ 'Tail' ].charAt(i) == 'C' ) $( '#read' + idx + 'pm' + i ).css({ 'color': colorC });
                            if( read[ 'Tail' ].charAt(i) == 'G' ) $( '#read' + idx + 'pm' + i ).css({ 'color': colorG });
                            if( read[ 'Tail' ].charAt(i) == 'T' ) $( '#read' + idx + 'pm' + i ).css({ 'color': colorT });

                            if( i != read[ 'Tail' ].length -1 )
                            {
                                $( '#read' + idx + 'pm' + i ).css({
                                    'width': spc_num + 'px',
                                    'float': 'left'
                                    });
                            }
                        }
                    });

                    $( '#chart' ).css({
                        'width': (( json[ 'Sequence' ].length + extra + lb_width )* spc_num ) + 'px',
                        'height': reads_count * read_height + seq_height + exp_height + 'px'
                        });

                    $( '#chart' ).append( \"<div id='line'></div>\" );

                    $( '#line' ).css({
                        'margin-top': fold_height + 'px',
                        'height': $( '#chart' ).height() + shift_top2 + 'px',
                        'width': '1px',
                        'position': 'absolute',
                        'top': '40px',
                        'left': '-1px',
                        'background-color': 'indianred',
                        'pointer-events': 'none'
                        });
                });

                var selectedfold = [];
                var selectedread = $( '#chart' );
                var selectedseed = $( '#chart' );

                $( document ).mousemove( function( event )
                {
                    if( event.pageY > shift_top1 )
                    {
                        $( '#line' ).css({
                            'left': ( event.pageX > $( '#chart' ).width() ? -1 : event.pageX ) + 'px'
                            });
                    }
                });

                $( function()
                {
                    $( document ).tooltip(
                    {
                        content: function()
                        {
                            var tip = '';
                            var data = Array();
                            selectedread = $( this );

                            if( selectedread.attr( 'title' ).substr( 0, 4 ) == 'fold' )
                                return tip;

                            if( selectedread.attr( 'title' ).substr( 0, 4 ) == 'expr' )
                            {
                                data = exprs[ selectedread.attr( 'title' ).substr( 4 ) ];
                                tip += 'PPM: ' + ( data[ 'PPM' ] * expr_max / 100 ).toFixed(2);
                                tip += '</br>Entropy: ' + ( data[ 'PPM' ] * expr_max / 100 ).toFixed(2);
                                tip += '</br>Entropy: ' + ( etpys[ selectedread.attr( 'title' ).substr( 4 ) ] ).toFixed(2);

                                if( data[ 'A' ] > 0.001 ) tip += '</br>Tail A: ' + ( data[ 'A' ] * 100 ).toFixed(2) + '%';
                                if( data[ 'C' ] > 0.001 ) tip += '</br>Tail C: ' + ( data[ 'C' ] * 100 ).toFixed(2) + '%';
                                if( data[ 'G' ] > 0.001 ) tip += '</br>Tail G: ' + ( data[ 'G' ] * 100 ).toFixed(2) + '%';
                                if( data[ 'T' ] > 0.001 ) tip += '</br>Tail T: ' + ( data[ 'T' ] * 100 ).toFixed(2) + '%';
                            }

                            if( selectedread.attr( 'title' ).substr( 0, 4 ) == 'Read' )
                            {
                                data = datas[ selectedread.attr( 'title' ).substr( 4 ) ];

                                var pos = '';
                                var MDseed = '';
                                var seq = sequ.substr( data[ 'Index' ], data[ 'Length' ]);
                                var MDarray = Array();
                                var TCarray = Array();

                                for( var i = 0; i < data[ 'MDtag' ].length; ++i )
                                {
                                    if( !$.isNumeric( data[ 'MDtag' ][i] ))
                                    {
                                        MDarray[ pos ] = data[ 'MDtag' ][i];
                                        pos = '';
                                    }
                                    else pos += data[ 'MDtag' ][i];
                                }

                                for( var i = 0; i < data[ 'TCtag' ].length; ++i )
                                {
                                    if( !$.isNumeric( data[ 'TCtag' ][i] ))
                                    {
                                        TCarray[ pos ] = data[ 'TCtag' ][i];
                                        pos = '';
                                    }
                                    else pos += data[ 'TCtag' ][i];
                                }

                                for( var i = 1; i < 8; ++i )
                                    if( i in MDarray ) MDseed += ( i -1 ) + MDarray[i];

                                $( '#selectbox' ).css({
                                    'left': ( data[ 'Index' ] * spc_num + 4 ) + 'px',
                                    'width': (( data[ 'Length' ] ) * spc_num ) + 'px',
                                    'display': 'block'
                                    });

                                selectedfold = [];

                                for( var i = 0; i < data[ 'Length' ]; ++i )
                                {
                                    selectedfold[i] = '#fold' + ( data[ 'Index' ] +i +1 );
                                    $( selectedfold[i] ).css({
                                        'background-color': 'Gold'
                                        });
                                }

                                tip += 'Annotation: ' + anno;
                                tip += ".( $Data_Array[0][5] != '.'? ( "'-' + data[ 'Arm' ]" ) : "''" ).";

                                tip += '_' + sequ.substr(( data[ 'Index' ] + 1 ), 7 );
                                tip += ( MDseed == '' ? '' : ( '|<b style=\"color:red;\">' + MDseed + '</b>' ));

                                tip += ( data[ 'isRMSK' ] == 'N' ? '' : ' (RMSK)' ) + '</br>';
                                tip += 'Sequence: ';

                                for( var i = 0; i < data[ 'Length' ]; ++i )
                                {
                                    if( i in MDarray )
                                    {
                                        tip += '<b style=\"color:red;\">' + MDarray[i] + '</b>';
                                    }
                                    else if( i in TCarray )
                                    {
                                        tip += '<b style=\"color:blue;\">' + TCarray[i] + '</b>';
                                    }
                                    else tip += seq[i];
                                }

                                tip += '</br>';

                                if( data[ 'MDtag' ] != '' ) tip += '<font color=red>MDtag: ' + data[ 'MDtag' ] + '</font></br>'; 
                                if( data[ 'TCtag' ] != '' ) tip += '<font color=blue>TCtag: ' + data[ 'TCtag' ] + '</font></br>'; 

                                tip += 'Length: ' + data[ 'Length' ];
                                tip += ( data[ 'Tail' ] == '' ? '</br>'
                                     : ( ' + ' + data[ 'Tail' ].length + '</br>Tail: ' + data[ 'Tail' ] + '</br>' ));

                                tip += 'PPM: ' + data[ 'PPM' ].toFixed(2); 
                            }

                            if( selectedread.attr( 'title' ).substr( 0, 4 ) == 'seqs' )
                            {
                                selectedseed = $( '#selectedseed' + ( selectedread.attr( 'title' ).substr( 4 ) - 5 ));

                                var idx = (( selectedread.attr( 'title' ).substr( 4 ) < 5 ) ? 0 :
                                          (( selectedread.attr( 'title' ).substr( 4 ) > sequ.length -3 ) ? ( sequ.length -7 ) :
                                           ( selectedread.attr( 'title' ).substr( 4 ) - 4 )));

                                $( '#selectbox' ).css({
                                    'left': ( idx + 0.225 ) * spc_num + 'px',
                                    'width': ( 7 * spc_num ) + 'px',
                                    'display': 'block'
                                    });

                                selectedseed.css({
                                    'border': '2px solid gold',
                                    'border-radius': '12px',
                                    });

                                if(( selectedread.attr( 'title' ).substr( 4 ) - 1 ) in seeds )
                                {
                                    tip += seeds[ selectedread.attr( 'title' ).substr( 4 ) - 1 ][ 'Name' ];
                                    tip += '</br>GMPM: ' + seeds[ selectedread.attr( 'title' ).substr( 4 ) - 1 ][ 'GMPM' ].toFixed(2);
                                    tip += '</br>GM: '   + seeds[ selectedread.attr( 'title' ).substr( 4 ) - 1 ][  'GM'  ].toFixed(2);
                                    tip += '</br>PM: '   + seeds[ selectedread.attr( 'title' ).substr( 4 ) - 1 ][  'PM'  ].toFixed(2);
                                }

                                return tip;
                            }

                            selectedread.css({
                                'border': '2px solid gold',
                                'border-radius': '12px'
                                });

                            return tip;
                        },

                        close: function( event, ui )
                        {
                            ui.tooltip.hover(
                                function ()
                                {
                                    $( this ).stop( true ).fadeTo( 400, 1 ); 

                                    if( selectedread.attr( 'title' ).substr( 0, 4 ) == 'Read' )
                                        $( '#selectbox' ).css({ 'display': 'block' });

                                    if( selectedread.attr( 'title' ).substr( 0, 4 ) == 'seqs' )
                                    {
                                        var idx = (( selectedread.attr( 'title' ).substr( 4 ) < 5 ) ? 0 :
                                                  (( selectedread.attr( 'title' ).substr( 4 ) > sequ.length -3 ) ? ( sequ.length -7 ) :
                                                   ( selectedread.attr( 'title' ).substr( 4 ) - 4 )));

                                        $( '#selectbox' ).css({
                                            'left': idx * spc_num + 'px',
                                            'width': 7 * spc_num + 'px',
                                            'display': 'block'
                                            });

                                        selectedseed.css({
                                            'border': '2px solid gold',
                                            'border-radius': '12px'
                                            });
                                    }
                                    else
                                    {
                                        selectedread.css({
                                            'border': '2px solid gold',
                                            'border-radius': '12px'
                                            });
                                    }
                                },
                                function ()
                                {
                                    $( this ).fadeOut( '400', function(){ $( this ).remove(); });

                                    if( selectedread.attr( 'title' ).substr( 0, 4 ) == 'Read' )
                                        $( '#selectbox' ).css({ 'display': 'none' });

                                    selectedread.css({ 'border': 'none' });
                                    selectedseed.css({ 'border': 'none' });
                                }
                            );

                            $( '#selectbox' ).css({ 'display': 'none' });
                            selectedread.css({ 'border': 'none' });
                            selectedseed.css({ 'border': 'none' });

                            for( var i = 0; i < selectedfold.length; ++i )
                            {
                                $( selectedfold[i] ).css({
                                    'background-color': 'transparent'
                                    });
                            }
                        }
                    });
                });

                </script>";
        }
    ?>
    </body>
</html>
