<? 
	Shell_Exec( 'rm /tmp/*' );
	$GMPM = $_POST['GMPM'];
	$FGMPM = $_POST['FGMPM'];
	$isLog = $_POST['isLog'];
	$RCPPM = $_POST['RCPPM'];
	$ForceY = $_POST['ForceY'];
	$Filter = $_POST['Filter'];
	$Bar_Pie = $_POST['Bar_Pie'];
	$MaxHight = $_POST['MaxHight'];
	$ForceMin = $_POST['ForceMin'];
	$ForceMax = $_POST['ForceMax'];
	$TSV_File = $_POST['TSV_File'];
	$TSV_File2 = $_POST['TSV_File2'];
	$Top_miRNA = $_POST['Top_miRNA'];
	$FilterMin = $_POST['FilterMin'];
	$FilterMax = $_POST['FilterMax'];
	$Chart_Type = $_POST['Chart_Type'];
	$isAbundant = $_POST['isAbundant'];
	$mirDistType = $_POST['mirDistType'];
	$miRNA_Select = $_POST['miRNA_Select'];
?>

<!DOCTYPE html>
<html>
	<meta charset='utf-8'>
	<body>

	<? 
#<!--================== Chart Type ====================-->

		echo '<script src=http://140.113.203.188/joye/ForAgoSorting/lib/d3.min.js></script>';
		echo '<link href=http://140.113.203.188/joye/ForAgoSorting/lib/svg0331.css rel=stylesheet type=text/css>';

		echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';

		echo '<select name=Chart_Type onchange=this.form.submit();>';
		echo '<option '; if($Chart_Type=='') echo 'selected'; echo '>Chart Type</option>';

		$Folder = Shell_Exec( 'ls -d */' );
		$Folder_List = Explode( "\n", $Folder );
		$Folder_Size = Count( $Folder_List );

		For( $i = 0; $i < $Folder_Size-1; ++$i )
		{
			echo '<option value='.$Folder_List[$i].' ';

			if( $Chart_Type == $Folder_List[$i] )
				echo 'selected ';

			echo '>'.$Folder_List[$i].'</option>';
		}
		echo "</select>
			</form>";

#<!--================== GMPM ====================-->

		if( $Chart_Type != 'ValPlot/' && $Chart_Type != 'LenPlus/' && $Chart_Type != 'MirTail/' )
		{
			echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';

			echo '<select name=GMPM onchange=this.form.submit();>';
			echo '<option '; if($GMPM=='') echo 'selected'; echo '>GM or PM</option>';

			$GMPM_List = array('GMPM', 'GM', 'PM', 'Tailing');
			$GMPM_Size = Count( $GMPM_List );

			For( $i = 0; $i < $GMPM_Size; ++$i )
			{
				echo '<option value='.$GMPM_List[$i].' ';

				if( $GMPM == $GMPM_List[$i] )
					echo 'selected ';

				echo '>' . $GMPM_List[$i] . '</option>';
			}
			echo "</select>
				<input type='hidden' name='FGMPM' value='$FGMPM' />
				<input type='hidden' name='isLog' value='$isLog' />
				<input type='hidden' name='Filter' value='$Filter' />
				<input type='hidden' name='Bar_Pie' value='$Bar_Pie' />
				<input type='hidden' name='MaxHight' value='$MaxHight' />
				<input type='hidden' name='ForceMin' value='$ForceMin' />
				<input type='hidden' name='ForceMax' value='$ForceMax' />
				<input type='hidden' name='TSV_File' value='$TSV_File' />
				<input type='hidden' name='TSV_File2' value='$TSV_File2' />
				<input type='hidden' name='Chart_Type' value='$Chart_Type' />
				<input type='hidden' name='isAbundant' value='$isAbundant' />
				<input type='hidden' name='mirDistType' value='$mirDistType' />
				</form>";

			if( $Chart_Type != 'MirDist/' && $Chart_Type != 'DotPlot/' )
			{
				echo '<script src=http://140.113.203.188/joye/ForAgoSorting/lib/nv.d3.min.js></script>';
				echo '<link href=http://140.113.203.188/joye/ForAgoSorting/lib/nv.d3.min.css rel=stylesheet type=text/css>';

#<!--================== Single TSV File ====================-->

				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';

				echo '<select name=TSV_File onchange=this.form.submit();>';
				echo '<option '; if($TSV_File=='') echo 'selected'; echo '>Select TSV</option>';

				$TSV = Shell_Exec( 'ls '.$Chart_Type.' | grep _'.$GMPM.'_ | grep .tsv' );
				$TSV_List = Explode( "\n", $TSV );
				$List_Size = Count( $TSV_List );

				For( $i = 0; $i < $List_Size-1; ++$i )
				{
					echo '<option value='.$TSV_List[$i].' ';

					if( $TSV_File == $TSV_List[$i] ) 
						echo 'selected ';

					echo '>'.$TSV_List[$i].'</option>';
				}
				echo "</select>
					<input type='hidden' name='Chart_Type' value='$Chart_Type' />
					<input type='hidden' name='GMPM' value='$GMPM' />
					<input type='hidden' name='Bar_Pie' value='$Bar_Pie' />
					</form>";

				if( $Chart_Type == 'BioType/' )
				{
					echo '<style>';
					echo 'text[ style="opacity: 0; text-anchor: middle;" ]{';
					echo 'opacity: 1 !important';
					echo '}';
					echo '</style>';

#<!--================== BioType ====================-->

					echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';

					echo '<input type=radio name=Bar_Pie value=Bar onchange=this.form.submit(); ';
					if( $Bar_Pie == 'Bar' ) echo 'checked';
					echo ' >Bar</input>';

					echo '<input type=radio name=Bar_Pie value=Pie onchange=this.form.submit(); ';
					if( $Bar_Pie == 'Pie' ) echo 'checked';
					echo ' >Pie</input>';

					echo "<input type='hidden' name='Chart_Type' value='$Chart_Type' />
						<input type='hidden' name='GMPM' value='$GMPM' />
						<input type='hidden' name='TSV_File' value='$TSV_File' />
						</form><br/>";

					$File = $Chart_Type.$TSV_File;
					$File_Tag = Explode( '_', $TSV_File );

					if( $Bar_Pie == 'Bar' )
					{
						if( $File_Tag[0] == 'All' )
						{
							echo "<svg id='bar'></svg>
								<script>
									d3.tsv( '$File', function( tsv_data ) {
										var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'BioType'; });
										var data = new Array();
							
										var id = 0;
										length.forEach( function(d) {
											var sample = new Object();
											var value = new Array();
							
											sample.key = String( [d] );
											sample.values = value;
										
											data[id] = sample;
											id++;
										});
							
										tsv_data.forEach( function(d) {
											id = 0;
							
											d.len = length.map( function(key) {
												var value = new Object();
												value.x = d['BioType'];
												value.y = Number(d[key]);
							
												data[id].values.push( value );
												id++;
											});
										});
							
										nv.addGraph({
											generate: function() {
												var width = nv.utils.windowSize().width,
													height = nv.utils.windowSize().height;
									
												var chart = nv.models.multiBarChart()
													.margin({'left':70})
													.staggerLabels(true)
													.width(width)
													.height(height);
									
												chart.dispatch.on('renderEnd', function(){
													console.log('Render Complete');
												});
								
												var svg = d3.select('#bar').datum(data);
												console.log('calling chart');
												svg.transition().duration(0).call(chart);
									
												return chart;
											},
											callback: function(graph) {
												nv.utils.windowResize(function() {
													var width = nv.utils.windowSize().width;
													var height = nv.utils.windowSize().height;
													graph.width(width).height(height);
									
													d3.select('#bar')
														.attr('width', width)
														.attr('height', height)
														.transition().duration(0)
														.call(graph);
									
												});
											}
										});
									});
								</script>";
						}
						else
						{
							echo "<svg id='bar'></svg>
								<script>
									d3.tsv( '$File', function( tsv_data ) {
										var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'BioType'; });
										var data = new Array();
							
										tsv_data.forEach( function(d) {
											d.len = length.map( function(key) {
												return {
													x: Number(key),
													y: +d[key]
												};
							
											});
											var sample = new Object();
											sample.key = d['BioType'];
											sample.values = d.len;
											data.push( sample );
										});
							
										nv.addGraph({
											generate: function() {
												var width = nv.utils.windowSize().width,
													height = nv.utils.windowSize().height;
									
												var chart = nv.models.multiBarChart()
													.margin({'left':70})
													.width(width)
													.height(height);
									
												chart.dispatch.on('renderEnd', function(){
													console.log('Render Complete');
												});
								
												var svg = d3.select('#bar').datum(data);
												console.log('calling chart');
												svg.transition().duration(0).call(chart);
									
												return chart;
											},
											callback: function(graph) {
												nv.utils.windowResize(function() {
													var width = nv.utils.windowSize().width;
													var height = nv.utils.windowSize().height;
													graph.width(width).height(height);
									
													d3.select('#bar')
														.attr('width', width)
														.attr('height', height)
														.transition().duration(0)
														.call(graph);
									
												});
											}
										});
									});
								</script>";
						}
					}

					if( $Bar_Pie == 'Pie' )
					{
						if( $File_Tag[0] == 'All' )
						{
							$inFile = File_get_contents( $File );
							$inFile_Lines = Explode( "\n", $inFile );
							$inFile_Line = Explode( "\t", $inFile_Lines[0] );
							$idCount = Count( $inFile_Line );

							For( $i = 1; $i < $idCount; ++$i )
							{
								echo '<svg id=pie'.$i.'></svg>';
							}

							echo "<script>
									d3.tsv( '$File', function( tsv_data ) {
										var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'BioType'; });
										var data = new Array();
							
										var id = 0;
										length.forEach( function(d) {
											var sample = new Object();
											var value = new Array();
							
											sample.key = String( [d] );
											sample.values = value;
							
											data[id] = sample;
											id++;
										});
							
										tsv_data.forEach( function(d) {
											id = 0;
							
											d.len = length.map( function(c) {
												var value = new Object();
												value.key = d['BioType'];
												value.y = d[c];
							
												data[id].values.push( value );
												id++;
											});
										});
							
										var svg_id = 0;
										data.forEach( function(d) {
											console.log( d.key );
											console.log( d.values );
							
											nv.addGraph({
												generate: function() {
													svg_id++;
										
													var chart = nv.models.pieChart()
														.x(function(c) { return c.key })
														.y(function(c) { return c.y })
														.donut(true)
														.width('500px')
														.height('500px')
														.cornerRadius(10)
														.showLabels(true);
								
													chart.title( d.key );
													chart.pie.donutLabelsOutside(true).donut(true);
								
													var svg = d3.select('#pie' + svg_id).datum( d.values );
													console.log('calling chart');
													svg.transition().duration(1200).call(chart);
								
													return chart;
												}
											});
							
										});
									});
								</script>";
						}
						else
						{
							echo '<br/>';

							$inFile = File_get_contents( $File );
							$inFile_Lines = Explode( "\n", $inFile );
							$idCount = Count( $inFile_Lines );

							For( $i = 1; $i < $idCount-1; ++$i )
							{
								echo '<svg id=pie'.$i.'></svg>';
							}

							echo "<script>
									d3.tsv( '$File', function( tsv_data ) {
										var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'BioType'; });
										var svg_id = 0;
							
										tsv_data.forEach( function(d) {
											d.len = length.map( function(key) {
												return {
													key: Number(key),
													y: +d[key]
												};
											});
							
											nv.addGraph({
												generate: function() {
													svg_id++;
							
													var chart = nv.models.pieChart()
														.x(function(c) { return c.key })
														.y(function(c) { return c.y })
														.donut(true)
														.width('500px')
														.height('500px')
														.cornerRadius(10)
														.showLabels(true);
										
													chart.title( d['BioType'] );
													chart.pie.donutLabelsOutside(true).donut(true);
								
													var svg = d3.select('#pie' + svg_id).datum( d.len );
													console.log('calling chart');
													svg.transition().duration(1200).call(chart);
								
													return chart;
												}
											});
										});
							
									});
							
								</script>";
						}
					}
				}

				if( $Chart_Type == 'LenDist/' )
				{

#<!--================== LenDist ====================-->

					echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';

					echo '<input type=radio name=Bar_Pie value=Bar onchange=this.form.submit(); ';
					if( $Bar_Pie == 'Bar' ) echo 'checked';
					echo ' >Bar</input>';

					echo '<input type=radio name=Bar_Pie value=Pie onchange=this.form.submit(); ';
					if( $Bar_Pie == 'Pie' ) echo 'checked';
					echo ' >Pie</input>';

					echo "
						<input type='hidden' name='Chart_Type' value='$Chart_Type' />
						<input type='hidden' name='GMPM' value='$GMPM' />
						<input type='hidden' name='TSV_File' value='$TSV_File' />
						</form><br/>";

					$File = $Chart_Type.$TSV_File;

					if( $Bar_Pie == 'Bar' )
					{
						echo "<svg id='bar'></svg>
							<script>
								d3.tsv( '$File', function( tsv_data ) {
								var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'miRNA'; });
								var data = new Array();
					
								tsv_data.forEach( function(d) {
									d.len = length.map( function(key) {
										return {
											x: Number(key),
											y: +d[key]
										};
					
									});
									var sample = new Object();
									sample.key = d['miRNA'];
									sample.values = d.len;
									data.push( sample );
								});
					
								nv.addGraph({
									generate: function() {
										var width = nv.utils.windowSize().width,
											height = nv.utils.windowSize().height;
							
										var chart = nv.models.multiBarChart()
											.margin({'left':70})
											.width(width)
											.height(height);
							
										chart.dispatch.on('renderEnd', function(){
											console.log('Render Complete');
										});
						
										var svg = d3.select('#bar').datum(data);
										console.log('calling chart');
										svg.transition().duration(0).call(chart);
							
										return chart;
									},
									callback: function(graph) {
										nv.utils.windowResize(function() {
											var width = nv.utils.windowSize().width;
											var height = nv.utils.windowSize().height;
											graph.width(width).height(height);
							
											d3.select('#bar')
												.attr('width', width)
												.attr('height', height)
												.transition().duration(0)
												.call(graph);
							
										});
									}
								});
							});
						</script>";
					}

					if( $Bar_Pie == 'Pie' )
					{
						$inFile = File_get_contents( $File );
						$inFile_Lines = Explode( "\n", $inFile );
						$inFile_Line = Explode( "\t", $inFile_Lines[0] );
						$idCount = Count( $inFile_Line );

						For( $i = 1; $i < $idCount; ++$i )
						{
							echo '<svg id=pie'.$i.'></svg>';
						}

						echo "<script>
								d3.tsv( '$File', function( tsv_data ) {
								var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'miRNA'; });
								var svg_id = 0;
					
								tsv_data.forEach( function(d) {
									d.len = length.map( function(key) {
										return {
											key: Number(key),
											y: +d[key]
										};
									});
					
									nv.addGraph({
										generate: function() {
											svg_id++;
								
											var chart = nv.models.pieChart()
												.x(function(c) { return c.key })
												.y(function(c) { return c.y })
												.donut(true)
												.width('500px')
												.height('500px')
												.cornerRadius(10)
												.showLabels(true);
								
											chart.title( d['miRNA'] );
											chart.pie.donutLabelsOutside(true).donut(true);
						
											var svg = d3.select('#pie' + svg_id).datum( d.len );
											console.log('calling chart');
											svg.transition().duration(1200).call(chart);
						
											return chart;
										}
									});
								});
					
							});
					
						</script>";
					}

				}
			}

			if( $Chart_Type == 'MirDist/' )
			{
				echo '<style>
						.x.axis path {
						display: none;
					}
					</style>';

#<!--================== Top miRNA Select ====================-->

				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';

				echo '<select name=Top_miRNA onchange=this.form.submit();>';
				echo '<option '; if($Top_miRNA=='') echo 'selected'; echo '>Select Top miRNA</option>';

				$TSV = Shell_Exec( 'ls '.$Chart_Type.' | grep _'.$GMPM.'_ | grep .tsv' );
				$TSV_List = Explode( "\n", $TSV );

				$Count_File = $Chart_Type . $TSV_List[0];
				$Count_inFile = File_get_contents( $Count_File );
				$Count_inFile_Lines = Explode( "\n", $Count_inFile );
				$Count_miRNA = Count( $Count_inFile_Lines )-2;
				$Count_miRNA10 = $Count_miRNA % 10;

				For( $i = 10; $i < $Count_miRNA; $i += 10 )
				{
					echo '<option value='.$i.' ';

					if( $i == $Top_miRNA )
						echo 'selected ';

					echo '>'.$i.'</option>';
				}

				echo '<option value='.$Count_miRNA.' ';

				if( $Count_miRNA == $Top_miRNA )
					echo 'selected ';

				echo '>'.$Count_miRNA.'</option>';

				echo "</select>";

#<!--================== mirDistType ====================-->

				echo '<select name=mirDistType onchange=this.form.submit();>';

				$mirDistType_List = array('ppm', '100%');
				$mirDistType_Size = Count( $mirDistType_List );

				For( $i = 0; $i < $mirDistType_Size; ++$i )
				{
					echo '<option value='.$mirDistType_List[$i].' ';

					if( $mirDistType == $mirDistType_List[$i] )
						echo 'selected ';

					echo '>' . $mirDistType_List[$i] . '</option>';
				}

				echo "</select>";

#<!--================== Max NU ====================-->

				if( $mirDistType == 'ppm' )
				{
					echo '<input type=text name=MaxHight size=8 value=';

					if( $MaxHight=='' )
						echo 'MaxHight';
					else
						echo $MaxHight;

					echo " onfocus=\"{this.value='';}\">";
				}

#<!--================== Multi TSV File ====================-->

				echo '<select name=TSV_File[] multiple=mutiple hight=>';

				$List_Size = Count( $TSV_List );

				For( $i = 0; $i < $List_Size-1; ++$i )
				{
					echo '<option value='.$TSV_List[$i].' ';

					foreach( $TSV_File as $TSVs )
						if( $TSVs == $TSV_List[$i] ) 
							echo 'selected ';

					echo '>'.$TSV_List[$i].'</option>';
				}
				echo "</select>
					<input type='hidden' name='Chart_Type' value='$Chart_Type' />
					<input type='hidden' name='GMPM' value='$GMPM' />
					<input type='submit' value='Submit' /> 
					</form><br/>";

#<!--================== MirDist ====================-->

				//========== svg view var set ==========

				$width = $Top_miRNA * 180;

				echo "<script>
						var margin = {top: 20, right: 60, bottom: 35, left: 40},
							width = $width - margin.left - margin.right, //180(per one miRNA)
							height = 300 - margin.top - margin.bottom;
						
						var x0 = d3.scale.ordinal()
							.rangeRoundBands([0, width], .3);
						
						var x1 = d3.scale.ordinal();
						
						var y = d3.scale.linear()
							.range([height, 0]);
						
						var color = d3.scale.ordinal()
							.range(['#F75C2F', '#E8B647', '#838A2D', '#66BAB7', '#6E75A4', '#72636E']);
						
						var xAxis = d3.svg.axis()
							.scale(x0)
							.orient('bottom');
						
						var yAxis = d3.svg.axis()
							.scale(y)
							.orient('left')
							.tickFormat(d3.format('.2s'));";

				//==================== svg ====================

				if( $mirDistType == 'ppm' )
				{
					if( $MaxHight == '' || $MaxHight == 'MaxHight' )
					{
						$MaxValue = 0;

						For( $i = 0; $i < Count( $TSV_File ); ++$i )
						{
							$File = $Chart_Type.$TSV_File[$i];
							$inFile = File_get_contents( $File );
							$inFile_Lines = Explode( "\n", $inFile );

							For( $j = 1; $j < Count( $inFile_Lines )-1; ++$j )
							{
								$inFile_Line = Explode( "\t", $inFile_Lines[$j] );

								For( $k = 1; $k < Count( $inFile_Line ); ++$k )
								{
									if( $inFile_Line[$k] >= $MaxValue )
										$MaxValue = $inFile_Line[$k];
								}
							}
						}
					}
					else
						$MaxValue = $MaxHight;
				}
				else
					$MaxValue = 100;

				For( $i = 0; $i < Count( $TSV_File ); ++$i )
				{
					$File = $Chart_Type.$TSV_File[$i];
					$inFile = File_get_contents( $File );
					$inFile_Lines = Explode( "\n", $inFile );

					$Temp = Tempnam( '/tmp', $TSV_File[$i] );
					$Ftemp = fopen( $Temp, 'w' );

					if( $mirDistType == 'ppm' )
					{
						For( $j = 0; $j < Count( $inFile_Lines )-1; ++$j )
						{
							if( $j >= $Top_miRNA+1 )
								break;

							fwrite( $Ftemp, $inFile_Lines[$j]."\n" );
						}
					}
					else
					{
						For( $j = 0; $j < Count( $inFile_Lines )-1; ++$j )
						{
							if( $j >= $Top_miRNA+1 )
								break;

							if( $j == 0 )
								fwrite( $Ftemp, $inFile_Lines[$j]."\n" );
							else
							{
								$inFile_Line_clm = Explode( "\t", $inFile_Lines[$j] );
								$Total_Value = 0;

								For( $k = 1; $k < Count( $inFile_Line_clm ); ++$k )
									$Total_Value += $inFile_Line_clm[ $k ];

								fwrite( $Ftemp, $inFile_Line_clm[0] );

								For( $k = 1; $k < Count( $inFile_Line_clm ); ++$k )
									fwrite( $Ftemp, "\t".Number_format( $inFile_Line_clm[ $k ]*100/$Total_Value, 0 ));

								fwrite( $Ftemp, "\n");
							}
						}
					}

					fclose( $Ftemp );

					echo "var svg$i = d3.select('body').append('svg')
							.attr('id', 'svg$i')
   	 					.attr('width', width + margin.left + margin.right)
   	 					.attr('height', height + margin.top + margin.bottom)
   	 					.append('g')
   	 					.attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');
   	 				
						d3.tsv('$Temp', function(error, data) {
   	 					var sample = d3.keys(data[0]).filter(function(key) { return key !== 'miRNA'; });
   	 				
   	 					data.forEach(function(d) {
   	 						 d.values = sample.map(function(name) { return {name: name, value: +d[name]}; });
   	 					});
   	 		
   	 					x0.domain(data.map(function(d) { return d.miRNA; }));
   	 					x1.domain(sample).rangeRoundBands([0, x0.rangeBand()]);
						y.domain([0, $MaxValue]);
   	 				
							svg$i.append('g')
   	 						.attr('class', 'x axis')
   	 						.attr('transform', 'translate(0,' + height + ')')
   	 						.call(xAxis);
   	 				
							svg$i.append('g')
   	 						.attr('class', 'y axis')
   	 						.call(yAxis)
   	 						.append('text')
   	 						.attr('transform', 'rotate(-90)')
   	 						.attr('y', 6)
   	 						.attr('dy', '.71em')
   	 						.style('text-anchor', 'end')
								.text('$TSV_File[$i]');
   	 				
							var barchat = svg$i.selectAll('.barchat')
   	 						.data(data)
   	 						.enter().append('g')
   	 						.attr('class', 'g')
   	 						.attr('transform', function(d) { return 'translate(' + x0(d.miRNA) + ',0)'; });
   	 		
   	 					var div = d3.select('body').append('div')
   	 						.attr('class', 'tooltip')
							.attr('id', 'Mir')
   	 						.style('opacity', 0);
   	 				
   	 					barchat.selectAll('rect')
   	 						.data(function(d) { return d.values; })
   	 						.enter().append('rect')
   	 						.attr('width', x1.rangeBand() - 5)
   	 						.attr('x', function(d) { return x1(d.name); })
   	 						.attr('y', function(d){
   	 							if( d.value < $MaxValue )
   	 								return y(d.value);
   	 							else
   	 								return y($MaxValue);
   	 						})
   	 						.attr('height', function(d){
   	 							if( d.value < $MaxValue )
   	 								return height - y(d.value);
   	 							else
   	 								return height - y($MaxValue);
   	 						})
   	 						.style('fill', function(d) { return color(d.name); })
   	 						.on('mouseover', function(d) {
   	 							div.transition()
   	 								.duration(200)
   	 								.style('opacity', .9);
   	 		
								div.html('<table align=center ><tr><th>Sample</th><th>Value</th></tr><tr><th>' +
										 d.name + '</th><th>' + d.value + '</th></tr>')
   	 								.style('left', (d3.event.pageX) + 'px')
   	 								.style('top', (d3.event.pageY - 28) + 'px');
   	 						})
   	 						.on('mouseout', function(d) {
   	 							div.transition()
   	 								.duration(500)
   	 								.style('opacity', 0);
   	 						});
   	 		
							var legend = svg$i.selectAll('.legend')
   	 						.data(sample.slice())
   	 						.enter().append('g')
   	 						.attr('class', 'legend')
   	 						.attr('transform', function(d, i) { return 'translate(0,' + i * 20 + ')'; });
   	 				
   	 					legend.append('rect')
   	 						.attr('x', width + 40)
   	 						.attr('width', 18)
   	 						.attr('height', 18)
   	 						.style('fill', color);
   	 				
   	 					legend.append('text')
   	 						.attr('x', width + 34)
   	 						.attr('y', 9)
   	 						.attr('dy', '.35em')
   	 						.style('text-anchor', 'end')
   	 						.text(function(d) { return d; });
   	 				});";
				}
				echo '</script>';
			}

			if( $Chart_Type == 'DotPlot/' )
			{

#<!--================== TSV File ====================-->
				
				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';

				$TSV = Shell_Exec( 'ls '.$Chart_Type.' | grep .tsv' );
				$TSV_List = Explode( "\n", $TSV );
				$List_Size = Count( $TSV_List );

				echo '<select name=TSV_File onchange=this.form.submit();>';
				echo '<option '; if($TSV_File=='') echo 'selected'; echo '>Select TSV</option>';

				For( $i = 0; $i < $List_Size-1; ++$i )
				{
					echo '<option value='.$TSV_List[$i].' ';

					if( $TSV_File == $TSV_List[$i] ) 
						echo 'selected ';

					echo '>'.$TSV_List[$i].'</option>';
				}
				echo '</select>';

				echo '<select name=TSV_File2 onchange=this.form.submit();>';
				echo '<option '; if($TSV_File2=='') echo 'selected'; echo '>Select TSV</option>';

				For( $i = 0; $i < $List_Size-1; ++$i )
				{
					echo '<option value='.$TSV_List[$i].' ';

					if( $TSV_File2 == $TSV_List[$i] ) 
						echo 'selected ';

					echo '>'.$TSV_List[$i].'</option>';
				}

				echo "</select>
					<input type='hidden' name='GMPM' value='$GMPM' />
					<input type='hidden' name='FGMPM' value='$FGMPM' />
					<input type='hidden' name='isLog' value='$isLog' />
					<input type='hidden' name='Filter' value='$Filter' />
					<input type='hidden' name='ForceMin' value='$ForceMin' />
					<input type='hidden' name='ForceMax' value='$ForceMax' />
					<input type='hidden' name='Chart_Type' value='$Chart_Type' />
					<input type='hidden' name='isAbundant' value='$isAbundant' />
					</form>";

#<!--================== is_Abundant ====================-->
				
				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';
				echo '<select name=isAbundant onchange=this.form.submit();>';
				
				$isAbundant_List = array('MostAbundant', 'AllmiRNA');
				$isAbundant_Size = Count( $isAbundant_List );
				if( $isAbundant == '' )
					$isAbundant = 'MostAbundant';
				
				For( $i = 0; $i < $isAbundant_Size; ++$i )
				{
					echo '<option value='.$isAbundant_List[$i].' ';
				
					if( $isAbundant == $isAbundant_List[$i] )
						echo 'selected ';
				
					echo '>' . $isAbundant_List[$i] . '</option>';
				}

				echo "</select>
					<input type='hidden' name='GMPM' value='$GMPM' />
					<input type='hidden' name='FGMPM' value='$FGMPM' />
					<input type='hidden' name='isLog' value='$isLog' />
					<input type='hidden' name='Filter' value='$Filter' />
					<input type='hidden' name='ForceMin' value='$ForceMin' />
					<input type='hidden' name='ForceMax' value='$ForceMax' />
					<input type='hidden' name='TSV_File' value='$TSV_File' />
					<input type='hidden' name='TSV_File2' value='$TSV_File2' />
					<input type='hidden' name='Chart_Type' value='$Chart_Type' />
					</form>";

#<!--================== isLog ====================-->

				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';

				echo '<select name=isLog onchange=this.form.submit();>';
				echo '<option '; if($isLog=='') echo 'selected'; echo 'value= >isLog</option>';

				$isLog_List = array(2, 4, 6, 8, 10);
				$isLog_Size = Count( $isLog_List );

				For( $i = 0; $i < $isLog_Size; ++$i )
				{
					echo '<option value='.$isLog_List[$i].' ';

					if( $isLog == $isLog_List[$i] )
						echo 'selected ';

					echo '>'.$isLog_List[$i].'</option>';
				}

				echo "</select>
					<input type='hidden' name='GMPM' value='$GMPM' />
					<input type='hidden' name='FGMPM' value='$FGMPM' />
					<input type='hidden' name='Filter' value='$Filter' />
					<input type='hidden' name='ForceMin' value='$ForceMin' />
					<input type='hidden' name='ForceMax' value='$ForceMax' />
					<input type='hidden' name='TSV_File' value='$TSV_File' />
					<input type='hidden' name='TSV_File2' value='$TSV_File2' />
					<input type='hidden' name='Chart_Type' value='$Chart_Type' />
					<input type='hidden' name='isAbundant' value='$isAbundant' />
					</form>";

#<!--================== ForceMin&Max ====================-->

				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';
				echo '<input type=text name=ForceMin size=5 value=';

				if( $ForceMin=='' )
					echo 'ForceMin';
				else
					echo $ForceMin;

				echo " onfocus=\"{this.value='';}\">";
				echo '<input type=text name=ForceMax size=5 value=';

				if( $ForceMax=='' )
					echo 'ForceMax';
				else
					echo $ForceMax;

				echo " onfocus=\"{this.value='';}\">";

#<!--================== Filter NU ====================-->

				echo '<input type=text name=Filter size=4 value=';

				if( $Filter=='' )
					echo 'Filter';
				else
					echo $Filter;

				echo " onfocus=\"{this.value='';}\">";

#<!--================== Filter GMPM ====================-->

				echo '<select name=FGMPM >';
				echo '<option '; if($FGMPM=='') echo 'selected'; echo ' value= >GM or PM</option>';

				$FGMPM_List = array('GMPM', 'GM', 'PM', 'Tailing');
				$FGMPM_Size = Count( $FGMPM_List );

				For( $i = 0; $i < $FGMPM_Size; ++$i )
				{
					echo '<option value='.$FGMPM_List[$i].' ';

					if( $FGMPM == $FGMPM_List[$i] )
						echo 'selected ';

					echo '>' . $FGMPM_List[$i] . '</option>';
				}

				echo "</select>
					<input type='hidden' name='GMPM' value='$GMPM' />
					<input type='hidden' name='isLog' value='$isLog' />
					<input type='hidden' name='Bar_Pie' value='$Bar_Pie' />
					<input type='hidden' name='TSV_File' value='$TSV_File' />
					<input type='hidden' name='TSV_File2' value='$TSV_File2' />
					<input type='hidden' name='Chart_Type' value='$Chart_Type' />
					<input type='hidden' name='isAbundant' value='$isAbundant' />
					<input type='submit' value='Submit' /> 
					</form><br/>";

#<!--================== Index & Header ====================-->

				$File1 = $Chart_Type.$TSV_File; 	//y
				$File2 = $Chart_Type.$TSV_File2;	//x
				$Files = array( $File1, $File2 );
				$Index = array();
				$Column = 0;
				$FColumn = 0;

				if( $isAbundant != 'MostAbundant' )
				{
					For( $l = 0; $l < Count( $Files ); ++$l )
					{
						$inFile = File_get_contents( $Files[$l] );
						$inFile_Lines = Explode( "\n", $inFile );

						For( $i = 0; $i < Count( $inFile_Lines )-1; ++$i )
						{
							$inFile_Line = Explode( "\t", $inFile_Lines[$i] );

							if( $i == 0 )
							{
								For( $j = 0; $j < Count( $inFile_Line ); ++$j )
								{
									if( $inFile_Line[$j] == $GMPM )
										$Column = $j;

									if( $inFile_Line[$j] == $FGMPM )
										$FColumn = $j;
								}
							}
							else
							{
								$miRNA_Seed = Explode( '*', $inFile_Line[0] );

								if( Count( $miRNA_Seed ) != 2 )
									Array_Push( $Index, $inFile_Line[0] );
								else
									Array_Push( $Index, $miRNA_Seed[1] );
							}
						}
					}
				}
				else
				{
					For( $l = 0; $l < Count( $Files ); ++$l )
					{
						$inFile = File_get_contents( $Files[$l] );
						$inFile_Lines = Explode( "\n", $inFile );

						For( $i = 0; $i < Count( $inFile_Lines )-1; ++$i )
						{
							$inFile_Line = Explode( "\t", $inFile_Lines[$i] );

							if( $i == 0 )
							{
								For( $j = 0; $j < Count( $inFile_Line ); ++$j )
								{
									if( $inFile_Line[$j] == $GMPM )
										$Column = $j;

									if( $inFile_Line[$j] == $FGMPM )
										$FColumn = $j;
								}

								Array_Push( $Index, 'miRNA' );
							}
							else
							{
								$miRNA_Seed = Explode( '*', $inFile_Line[0] );

								if( Count( $miRNA_Seed ) != 2 )
								{}
								else
									Array_Push( $Index, $miRNA_Seed[1] );
							}
						}
					}
				}

				Sort( $Index, SORT_STRING );
				$uIndex = Array_Unique( $Index, SORT_STRING );

#<!--================== Read File & Log ====================-->

				$Anno_Value1 = array();
				$Anno_Value2 = array();
				$Filter_Array1 = array();
				$Filter_Array2 = array();
				$Sample_Name = array();
				$FSample_Name = array();

				For( $l = 0; $l < Count( $Files ); ++$l )
				{
					$inFile = File_get_contents( $Files[$l] );
					$inFile_Lines = Explode( "\n", $inFile );

					For( $i = 0; $i < Count( $inFile_Lines )-1; ++$i )
					{
						$inFile_Line = Explode( "\t", $inFile_Lines[$i] );

						if( $i == 0 )
						{
							if( $l == 0 )
							{
								$Anno_Value1['miRNA'] = $inFile_Line[0];
								$Filter_Array1['miRNA']=$inFile_Line[0].'_F'.$FGMPM;
							}
							else
							{
								$Anno_Value2['miRNA'] = $inFile_Line[0];
								$Filter_Array2['miRNA']=$inFile_Line[0].'_F'.$FGMPM;
							}

							Array_Push( $Sample_Name, $inFile_Line[0] );
							Array_Push( $FSample_Name, $inFile_Line[0].'_F'.$FGMPM );
						}
						else
						{
							$miRNA_Seed = Explode( '*', $inFile_Line[0] );
							$Value = $inFile_Line[$Column];
							$FValue = $inFile_Line[$FColumn];

							if( $isLog != '' && $Value != 0 )
								$Value = ( Log($inFile_Line[$Column]) / Log($isLog) );

							if( $isLog != '' && $FValue != 0 )
								$FValue = ( Log($inFile_Line[$FColumn]) / Log($isLog) );

							if( Count( $miRNA_Seed ) != 2 )
							{
								if( $l == 0 )
								{
									$Anno_Value1[$inFile_Line[0]] = Round($Value,2);
									$Filter_Array1[$inFile_Line[0]]=Round($FValue,2);
								}
								else
								{
									$Anno_Value2[$inFile_Line[0]] = Round($Value,2);
									$Filter_Array2[$inFile_Line[0]]=Round($FValue,2);
								}
							}
							else
							{
								if( $l == 0 )
								{
									$Anno_Value1[$miRNA_Seed[1]] = Round($Value,2);
									$Filter_Array1[$miRNA_Seed[1]]=Round($FValue,2);
								}
								else
								{
									$Anno_Value2[$miRNA_Seed[1]] = Round($Value,2);
									$Filter_Array2[$miRNA_Seed[1]]=Round($FValue,2);
								}
							}
						}
					}
				}

#<!--================== Filter & Temp ====================-->

				$Temp = Tempnam( '/tmp', $GMPM.'_'.$Sample_Name[0].'_'.$Sample_Name[1].'_'.$isAbundant.'_'.$isLog.'_'.$Filter.'_'.$FGMPM );
				$Ftemp = Fopen( $Temp, 'w' );
				$MaxAxis = 0;

				Fwrite( $Ftemp, 'miRNA'."\t".
						$Anno_Value1['miRNA']."\t".
						$Anno_Value2['miRNA']."\t".
						$Filter_Array1['miRNA']."\t".
						$Filter_Array2['miRNA']."\n" );

				For( $i = 0; $i < Count( $Index ); ++$i )
				{
					if( $uIndex[$i] != '' && $uIndex[$i] != 'miRNA' )
					{
						if( $Anno_Value1[$uIndex[$i]] == '' )
						{
							$Anno_Value1[$uIndex[$i]]  = 0;
							$Filter_Array1[$uIndex[$i]]= 0;
						}

						if( $Anno_Value2[$uIndex[$i]] == '' )
						{
							$Anno_Value2[$uIndex[$i]]  = 0;
							$Filter_Array2[$uIndex[$i]]= 0;
						}

						if( $Filter != 'Filter' && $FGMPM != '' )
							if( $Filter_Array1[$uIndex[$i]] < $Filter && $Filter_Array2[$uIndex[$i]] < $Filter )
								continue;

						Fwrite( $Ftemp, $uIndex[$i]."\t".
										$Anno_Value1[$uIndex[$i]]."\t".
										$Anno_Value2[$uIndex[$i]]."\t".
										$Filter_Array1[$uIndex[$i]]."\t".
										$Filter_Array2[$uIndex[$i]]."\n" );

						if( $MaxAxis < $Anno_Value1[$uIndex[$i]] )
							$MaxAxis = $Anno_Value1[$uIndex[$i]];

						if( $MaxAxis < $Anno_Value2[$uIndex[$i]] )
							$MaxAxis = $Anno_Value2[$uIndex[$i]];
					}
				}

				Fclose( $Ftemp );

				if( $FGMPM == '' )
					$FGMPM = 'Filter';

				if( $ForceMin == 'ForceMin' || $ForceMin == '' )
					$ForceMin = 0;

				if( $ForceMax == 'ForceMax' || $ForceMax == '' )
					$ForceMax = $MaxAxis;

#<!--================== DotPlot ====================-->
			
				echo "<script>
					var svg_width  = window.innerWidth;
					var svg_height = window.innerHeight;
			
					var margin = {top: 20, right: 60, bottom: 60, left: 60}
						width = svg_width - margin.left - margin.right,
						height = svg_height - margin.top - margin.bottom;
			
					var x = d3.scale.linear()
						.range([0, width]);
			
					var y = d3.scale.linear()
						.range([height, 0]);
			
					var color = d3.scale.category10();
			
					var xAxis = d3.svg.axis()
						.scale(x)
						.orient('bottom');
			
					var yAxis = d3.svg.axis()
						.scale(y)
						.orient('left');
			
					var svg = d3.select('body').append('svg')
						.attr('width', width)
						.attr('height', height)
					.append('g')
						.attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');
			
					d3.tsv('$Temp', function(error, data) {
						data.forEach(function(d) {
							d.$Sample_Name[1] = +d.$Sample_Name[1];
							d.$Sample_Name[0] = +d.$Sample_Name[0];
						});
			
						var xExtent = d3.extent(data, function(d) { return d.$Sample_Name[1]; });
						var yExtent = d3.extent(data, function(d) { return d.$Sample_Name[0]; });
			
						xExtent[0] = $ForceMin;
						yExtent[0] = $ForceMin;
			
						xExtent[1] = $ForceMax;
						yExtent[1] = $ForceMax;
			
						x.domain( xExtent ).nice();
						y.domain( yExtent ).nice();
			
						svg.append('g')
							.attr('class', 'x axis')
							.attr('transform', 'translate(0,' + height + ')')
							.call(xAxis)
							.append('text')
							.attr('class', 'label')
							.attr('x', width)
							.attr('y', -6)
							.style('text-anchor', 'end')
							.text('$Sample_Name[1]');
			
						svg.append('g')
							.attr('class', 'y axis')
							.call(yAxis)
							.append('text')
							.attr('class', 'label')
							.attr('transform', 'rotate(-90)')
							.attr('y', 6)
							.attr('dy', '.71em')
							.style('text-anchor', 'end')
							.text('$Sample_Name[0]');
			
						svg.append('g')
							.attr('class', 'xyline')
							.append('line')
							.attr('x1', x($ForceMin) )
							.attr('y1', y($ForceMin) )
							.attr('x2', x($ForceMax) )
							.attr('y2', y($ForceMax) )
							.style('stroke','rgb(255,0,0)')
							.style('stroke-width','1');
			
						var div = d3.select('body').append('div')
							.attr('class', 'tooltip')
							.attr('id', 'Dot')
							.style('opacity', 0);
			
						svg.selectAll('.dot')
							.data(data)
							.enter().append('circle')
							.attr('class', 'dot')
							.attr('r', 5)
							.attr('cx', function(d) { return x(d.$Sample_Name[1]); })
							.attr('cy', function(d) { return y(d.$Sample_Name[0]); })
							.style('fill', function(d) { return color(d.species); })
							.on('mouseover', function(d) {
								div.transition()
									.duration(200)
									.style('opacity', .9);
			
								div.html('<table align=center ><tr><th colspan=3 >' + d.miRNA +
										 '</th></tr><tr><th>Sample</th><th>$GMPM</th><th>$FGMPM</th></tr><tr><th>$Sample_Name[0]</th><th>' +
										 d.$Sample_Name[0] + '</th><th>' + d.$Sample_Name[0]_F$FGMPM + '</th></tr><tr><th>$Sample_Name[1]</th><th>' +
										 d.$Sample_Name[1] + '</th><th>' + d.$Sample_Name[1]_F$FGMPM + '</th></tr>')
									.style('left', (d3.event.pageX) + 'px')
									.style('top', (d3.event.pageY - 28) + 'px');
							})
							.on('mouseout', function(d) {
								div.transition()
									.duration(500)
									.style('opacity', 0);
							});
			
							window.onresize = function(){
							};
					});
				</script>";
			}
		}
		else
		{
			if( $Chart_Type == 'ValPlot/' )
			{
				echo '	<script src=http://140.113.203.188/joye/ForAgoSorting/lib/d3.min.js ></script>
						<script src=http://140.113.203.188/joye/ForAgoSorting/lib/slickgrid/jquery-1.7.min.js ></script>
						<script src=http://140.113.203.188/joye/ForAgoSorting/lib/slickgrid/jquery.event.drag-2.0.min.js ></script>
						<script src=http://140.113.203.188/joye/ForAgoSorting/lib/slickgrid/slick.core.js ></script>
						<script src=http://140.113.203.188/joye/ForAgoSorting/lib/slickgrid/slick.grid.js ></script>
						<script src=http://140.113.203.188/joye/ForAgoSorting/lib/slickgrid/slick.dataview.js ></script>
						<script src=http://140.113.203.188/joye/ForAgoSorting/lib/d3.parcoords.js ></script>
						<script src=http://140.113.203.188/joye/ForAgoSorting/lib/divgrid.js ></script>
						<link rel=stylesheet type=text/css href=http://140.113.203.188/joye/ForAgoSorting/lib/slickgrid/slick.grid.css />
						<link rel=stylesheet type=text/css href=http://140.113.203.188/joye/ForAgoSorting/lib/slickgrid/jquery-ui-1.8.16.custom.css />
						<link rel=stylesheet type=text/css href=http://140.113.203.188/joye/ForAgoSorting/lib/slickgrid/examples.css />
						<link rel=stylesheet type=text/css href=http://140.113.203.188/joye/ForAgoSorting/lib/d3.parcoords.css />
						<link rel=stylesheet type=text/css href=http://140.113.203.188/joye/ForAgoSorting/lib/style.css />';

#<!--================== TSV File ====================-->

				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';

				echo '<select name=TSV_File onchange=this.form.submit();>';
				echo '<option '; if($TSV_File=='') echo 'selected'; echo '>Select TSV</option>';

				$TSV = Shell_Exec( 'ls '.$Chart_Type.' | grep .tsv' );
				$TSV_List = Explode( "\n", $TSV );
				$List_Size = Count( $TSV_List );

				For( $i = 0; $i < $List_Size-1; ++$i )
				{
					echo '<option value='.$TSV_List[$i].' ';

					if( $TSV_File == $TSV_List[$i] ) 
						echo 'selected ';

					echo '>'.$TSV_List[$i].'</option>';
				}
				echo "</select>
					<input type='hidden' name='FilterMin' value='$FilterMin' />
					<input type='hidden' name='FilterMax' value='$FilterMax' />
					<input type='hidden' name='Chart_Type' value='$Chart_Type' />
					</form>";

#<!--================== is_Abundant ====================-->
				
				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';
				echo '<select name=isAbundant onchange=this.form.submit();>';
				
				$isAbundant_List = array('MostAbundant', 'AllmiRNA');
				$isAbundant_Size = Count( $isAbundant_List );
				if( $isAbundant == '' )
					$isAbundant = 'MostAbundant';
				
				For( $i = 0; $i < $isAbundant_Size; ++$i )
				{
					echo '<option value='.$isAbundant_List[$i].' ';
				
					if( $isAbundant == $isAbundant_List[$i] )
						echo 'selected ';
				
					echo '>' . $isAbundant_List[$i] . '</option>';
				}

				echo "</select>
					<input type='hidden' name='FGMPM' value='$FGMPM' />
					<input type='hidden' name='isLog' value='$isLog' />
					<input type='hidden' name='Filter' value='$Filter' />
					<input type='hidden' name='TSV_File' value='$TSV_File' />
					<input type='hidden' name='FilterMin' value='$FilterMin' />
					<input type='hidden' name='FilterMax' value='$FilterMax' />
					<input type='hidden' name='Chart_Type' value='$Chart_Type' />
					</form>";

#<!--================== isLog ====================-->

				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';

				echo '<select name=isLog onchange=this.form.submit();>';
				echo '<option '; if($isLog=='') echo 'selected'; echo 'value= >isLog</option>';

				$isLog_List = array(2, 4, 6, 8, 10);
				$isLog_Size = Count( $isLog_List );

				For( $i = 0; $i < $isLog_Size; ++$i )
				{
					echo '<option value='.$isLog_List[$i].' ';

					if( $isLog == $isLog_List[$i] )
						echo 'selected ';

					echo '>'.$isLog_List[$i].'</option>';
				}

				echo "</select>
					<input type='hidden' name='FGMPM' value='$FGMPM' />
					<input type='hidden' name='Filter' value='$Filter' />
					<input type='hidden' name='TSV_File' value='$TSV_File' />
					<input type='hidden' name='FilterMin' value='$FilterMin' />
					<input type='hidden' name='FilterMax' value='$FilterMax' />
					<input type='hidden' name='Chart_Type' value='$Chart_Type' />
					<input type='hidden' name='isAbundant' value='$isAbundant' />
					</form>";

#<!--================== Filter Min & Max ====================-->

				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';
				echo '<input type=text name=FilterMin size=5 value=';

				if( $FilterMin=='' )
					echo 'FilterMin';
				else
					echo $FilterMin;

				echo " onfocus=\"{this.value='';}\">";
				echo '<input type=text name=FilterMax size=5 value=';

				if( $FilterMax=='' )
					echo 'FilterMax';
				else
					echo $FilterMax;

				echo " onfocus=\"{this.value='';}\">";

#<!--================== Filter NU ====================-->

				echo '<input type=text name=Filter size=4 value=';

				if( $Filter=='' )
					echo 'Filter';
				else
					echo $Filter;

				echo " onfocus=\"{this.value='';}\">";

#<!--================== Filter GMPM ====================-->

				echo '<select name=FGMPM >';
				echo '<option '; if($FGMPM=='') echo 'selected'; echo ' value= >GM or PM</option>';

				For( $i = 0; $i < $List_Size-1; ++$i )
				{
					echo '<option value='.$TSV_List[$i].' ';

					if( $FGMPM == $TSV_List[$i] ) 
						echo 'selected ';

					echo '>'.$TSV_List[$i].'</option>';
				}

				echo "
					<input type='hidden' name='isLog' value='$isLog' />
					<input type='hidden' name='TSV_File' value='$TSV_File' />
					<input type='hidden' name='Chart_Type' value='$Chart_Type' />
					<input type='hidden' name='isAbundant' value='$isAbundant' />
					<input type='submit' value='Submit' /> 
					</form><br/>";

#<!--================== Read File ====================-->

				if( $TSV_File != '' )
				{
					if( $FGMPM != '' && $Filter != '' )
					{
						$FFile = $Chart_Type.$FGMPM;
						$FinFile = File_get_contents( $FFile );
						$FinFile_Lines = Explode( "\n", $FinFile );
						$FArray = Array();
					}

					$File = $Chart_Type.$TSV_File;
					$inFile = File_get_contents( $File );
					$inFile_Lines = Explode( "\n", $inFile );

					$Temp = Tempnam( '/tmp', $TSV_File );
					$Ftemp = Fopen( $Temp, 'w' );

					$Sample_Nu = Count( Explode( "\t", $inFile_Lines[0] ))-1;
					$Value_Max = 0;
					$Value_Min = 0;

					For( $i = 0; $i < Count( $inFile_Lines )-1; ++$i )
					{
						if( $FGMPM != '' && $Filter != '' && $i != 0 )
						{
							$Check_Filter = Array();
							$FinFile_Line = Explode( "\t", $FinFile_Lines[$i] );

							For( $j = 1; $j < Count( $FinFile_Line ); ++$j )
							{
								$FValue = $FinFile_Line[$j];

								if( $isLog != '' && $FValue != 0 )
									$FValue = ( Log($FinFile_Line[$j]) / Log($isLog) );

								if( $FValue <= $Filter )
									Array_Push( $Check_Filter, $FValue );
							}

							if( Count( $Check_Filter ) == Count( $FinFile_Line )-1 )
								Continue;
						}


						$Check_Filter = 'keep';

						if( $i != 0 )
						{
							$Check_Min = Array();
							$Check_Max = Array();

							$inFile_Line = Explode( "\t", $inFile_Lines[$i] );

							For( $j = 1; $j < Count( $inFile_Line ); ++$j )
							{
								$Value = $inFile_Line[$j];

								if( $isLog != '' && $Value != 0 )
									$Value = ( Log($inFile_Line[$j]) / Log($isLog) );

								if( $FilterMin != 'FilterMin' && $FilterMin != '' )
									if( $Value <= $FilterMin )
										Array_Push( $Check_Min, $Value );

								if( $FilterMax != 'FilterMax' && $FilterMax != '' )
									if( $Value >= $FilterMax )
										Array_Push( $Check_Max, $Value );
							}

							if( Count( $Check_Min ) == Count( $inFile_Line )-1 || Count( $Check_Max ) == Count( $inFile_Line )-1 )
								$Check_Filter = 'notkeep';


							if( $isAbundant == 'MostAbundant' )
							{
								$miRNA = Explode( '*', $inFile_Line[0] );

								if( Count( $miRNA ) != 2 )
									$Check_Filter = 'notkeep';
							}

							if( $Check_Filter == 'keep' )
							{
								For( $j = 1; $j < Count( $inFile_Line ); ++$j )
								{
									$Value = $inFile_Line[$j];

									if( $isLog != '' && $Value != 0 )
										$Value = ( Log($inFile_Line[$j]) / Log($isLog) );

									if( $Value_Max <= $Value )
										$Value_Max =  $Value;

									if( $Value_Min >= $Value )
										$Value_Min =  $Value;
								}
							}
						}

						if( $Check_Filter == 'keep' )
						{
							if( $isLog != '' && $i != 0)
							{
								$inFile_Line = Explode( "\t", $inFile_Lines[$i] );
								Fwrite( $Ftemp, $inFile_Line[0] );

								For( $j = 1; $j < Count( $inFile_Line ); ++$j )
								{
									if( $inFile_Line[$j] != 0 )
									{	
										$Value = ( Log($inFile_Line[$j]) / Log($isLog) );
										Fwrite( $Ftemp, "\t".$Value );
									}
									else
										Fwrite( $Ftemp, "\t".'0' );
								}

								Fwrite( $Ftemp, "\n" );
							}
							else
								Fwrite( $Ftemp, $inFile_Lines[$i]."\n" );
						}
					}

					Fwrite( $Ftemp, 'MinValue' );
					For( $i = 0; $i < $Sample_Nu; ++$i )
						Fwrite( $Ftemp, "\t".$Value_Min );
					Fwrite( $Ftemp, "\n" );

					Fwrite( $Ftemp, 'MaxValue' );
					For( $i = 0; $i < $Sample_Nu; ++$i )
						Fwrite( $Ftemp, "\t".$Value_Max );
					Fwrite( $Ftemp, "\n" );

					fclose( $Ftemp );

#<!--================== ValPlot ====================-->
	
					echo "<div id='example' class='parcoords' style='height:240px;'></div>
						<div id='grid'></div>
						<script id='brushing'>
							var parcoords = d3.parcoords()('#example')
									.alpha(0.4)
									.mode('queue') // progressive rendering
									.height(d3.max([document.body.clientHeight-326, 220]))
									.margin({
										top: 36,
										left: 0,
										right: 0,
										bottom: 16
									});
							
							// load csv file and create the chart
							d3.tsv('$Temp', function(data) {
								// slickgrid needs each data element to have an id
								data.forEach(function(d,i) { d.id = d.id || i; });
							
								parcoords
									.data(data)
									.hideAxis(['miRNA'])
									.hideAxis(['id'])
									.render()
									.reorderable()
									.brushMode('1D-axes');
							
								// setting up grid
								var column_keys = d3.keys(data[0]);
								var columns = column_keys.map(function(key,i) {
									return {
										id: key,
										name: key,
										field: key,
										sortable: true
									}
								});
							
								var options = {
									enableCellNavigation: true,
									enableColumnReorder: false,
									multiColumnSort: false
								};
							
								var dataView = new Slick.Data.DataView();
								var grid = new Slick.Grid('#grid', dataView, columns, options);
							
								// wire up model events to drive the grid
								dataView.onRowCountChanged.subscribe(function (e, args) {
									grid.updateRowCount();
									grid.render();
								});
							
								dataView.onRowsChanged.subscribe(function (e, args) {
									grid.invalidateRows(args.rows);
									grid.render();
								});
							
								// column sorting
								var sortcol = column_keys[0];
								var sortdir = 1;
							
								function comparer(a, b) {
									var x = a[sortcol], y = b[sortcol];
									return (x == y ? 0 : (x > y ? 1 : -1));
								}
								
								// click header to sort grid column
								grid.onSort.subscribe(function (e, args) {
									sortdir = args.sortAsc ? 1 : -1;
									sortcol = args.sortCol.field;
							
									if ($.browser.msie && $.browser.version <= 8) {
										dataView.fastSort(sortcol, args.sortAsc);
									} else {
										dataView.sort(comparer, args.sortAsc);
									}
								});
							
								// highlight row in chart
								grid.onMouseEnter.subscribe(function(e,args) {
									var i = grid.getCellFromEvent(e).row;
									var d = parcoords.brushed() || data;
									parcoords.highlight([d[i]]);
								});
								grid.onMouseLeave.subscribe(function(e,args) {
									parcoords.unhighlight();
								});
							
								// fill grid with data
								gridUpdate(data);
							
								// update grid on brush
								parcoords.on('brush', function(d) {
									gridUpdate(d);
								});
							
								function gridUpdate(data) {
									dataView.beginUpdate();
									dataView.setItems(data);
									dataView.endUpdate();
								};
							});
						</script>";
				}
			}

			if( $Chart_Type == 'LenPlus/' )
			{
				echo '<script src=http://140.113.203.188/joye/ForAgoSorting/lib/nv.d3.min.js></script>';
				echo '<link href=http://140.113.203.188/joye/ForAgoSorting/lib/nv.d3.min.css rel=stylesheet type=text/css>';
				echo '<script src=http://140.113.203.188/joye/ForAgoSorting/lib/head.min.js ></script>';
				echo '<script src=http://140.113.203.188/joye/ForAgoSorting/lib/reveal.min.js></script>';
				echo '<link rel=stylesheet href=http://140.113.203.188/joye/ForAgoSorting/lib/reveal.min.css>';

#<!--================== ReadCount & PPM ====================-->

				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';
				echo '<select name=RCPPM onchange=this.form.submit();>';

				$RCPPM_List = array('ppm', 'readcount');
				$RCPPM_Size = Count( $RCPPM_List );
				if( $RCPPM == '' )
					$RCPPM = 'ppm';

				For( $i = 0; $i < $RCPPM_Size; ++$i )
				{
					echo '<option value='.$RCPPM_List[$i].' ';

					if( $RCPPM == $RCPPM_List[$i] )
						echo 'selected ';

					echo '>' . $RCPPM_List[$i] . '</option>';
				}

				echo "</select>
					<input type='hidden' name='ForceY' value='$ForceY' />
					<input type='hidden' name='TSV_File' value='$TSV_File' />
					<input type='hidden' name='Chart_Type' value='$Chart_Type' />
					</form>";

#<!--================== ForceY ====================-->

				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';
				echo '<input type=text name=ForceY size=3 value=';

				if( $ForceY=='' )
					echo 'Hight';
				else
					echo $ForceY;

				echo " onfocus=\"{this.value='';}\">";

				echo "</select>
					<input type='hidden' name='RCPPM' value='$RCPPM' />
					<input type='hidden' name='TSV_File' value='$TSV_File' />
					<input type='hidden' name='Chart_Type' value='$Chart_Type' />
					<input type='submit' value='Submit' /> 
					</form>";

#<!--================== TSV File ====================-->

				$TSV = Shell_Exec( 'ls '.$Chart_Type.' | grep '.$RCPPM.'.tsv' );
				$TSV_List = Explode( "\n", $TSV );
				$List_Size = Count( $TSV_List );

				echo '<div class=reveal><div class=slides>';

				For( $i = 0; $i < $List_Size-1; ++$i )
				{
					$AGO = Explode( '_', $TSV_List[$i]);
					echo '<section>';
					echo '<h2>'.$AGO[0].'</h2>';
					echo '<svg id=bar'.$i.'></svg>';
					echo '</section>';
				}

#<!--================== Slide Set ====================-->

				echo '</div></div>';
				echo "<script>
						Reveal.initialize({
							center: false,
							transition: Reveal.getQueryHash().transition || 'fade', // default/cube/page/concave/zoom/linear/fade/none
						});
					</script>";

#<!--================== Len Chart ====================-->

				For( $i = 0; $i < $List_Size-1; ++$i )
				{
					$File = $Chart_Type.$TSV_List[$i];

					echo "<script>
							d3.tsv( '$File', function( tsv_data ) {
							var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'Length'; });
							var data = new Array();
				
							tsv_data.forEach( function(d) {
								d.len = length.map( function(key) {
									return {
										x: Number(key),
										y: +d[key]
									};
				
								});
								var sample = new Object();
								sample.key = d['Length'];
								sample.values = d.len;
								data.push( sample );
							});
				
							nv.addGraph({
								generate: function() {
									var width = nv.utils.windowSize().width,
										height = nv.utils.windowSize().height;
						
									var chart = nv.models.multiBarChart()";

						if( $ForceY != '' && $ForceY != 'Hight' )
							echo "		.forceY([$ForceY,0])";

					echo "				.width(width)
										.height(height)
										.color( ['#000000', '#FF0000', '#0000FF', '#FFBF00', '#088A08', '#6E6E6E'])
										.stacked(true);
					
									var svg = d3.select('#bar$i').datum(data);
									console.log('calling chart');
									svg.transition().duration(0).call(chart);
						
									return chart;
								}
							});
						});
					</script>";

				}
			}

			if( $Chart_Type == 'MirTail/' )
			{
				echo '<script src=http://140.113.203.188/joye/ForAgoSorting/lib/nv.d3.min.js></script>';
				echo '<link href=http://140.113.203.188/joye/ForAgoSorting/lib/nv.d3.min.css rel=stylesheet type=text/css>';

#<!--================== TSV File ====================-->

				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';

				echo '<select name=TSV_File onchange=this.form.submit();>';
				echo '<option '; if($TSV_File=='') echo 'selected'; echo '>Select TSV</option>';

				$TSV = Shell_Exec( 'ls '.$Chart_Type.' | grep .tsv' );
				$TSV_List = Explode( "\n", $TSV );
				$List_Size = Count( $TSV_List );

				For( $i = 0; $i < $List_Size-1; ++$i )
				{
					echo '<option value='.$TSV_List[$i].' ';

					if( $TSV_File == $TSV_List[$i] ) 
						echo 'selected ';

					echo '>'.$TSV_List[$i].'</option>';
				}
				echo "</select>
					<input type='hidden' name='ForceY' value='$ForceY' />
					<input type='hidden' name='Chart_Type' value='$Chart_Type' />
					<input type='hidden' name='miRNA_Select' value='$miRNA_Select' />
					</form>";

#<!--================== is_Abundant ====================-->
				
				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';
				echo '<select name=isAbundant onchange=this.form.submit();>';
				
				$isAbundant_List = array('MostAbundant', 'AllmiRNA');
				$isAbundant_Size = Count( $isAbundant_List );
				if( $isAbundant == '' )
					$isAbundant = 'MostAbundant';
				
				For( $i = 0; $i < $isAbundant_Size; ++$i )
				{
					echo '<option value='.$isAbundant_List[$i].' ';
				
					if( $isAbundant == $isAbundant_List[$i] )
						echo 'selected ';
				
					echo '>' . $isAbundant_List[$i] . '</option>';
				}

				echo "</select>
					<input type='hidden' name='ForceY' value='$ForceY' />
					<input type='hidden' name='TSV_File' value='$TSV_File' />
					<input type='hidden' name='Chart_Type' value='$Chart_Type' />
					<input type='hidden' name='miRNA_Select' value='$miRNA_Select' />
					</form>";

#<!--================== miRNA Select ====================-->

				$File = $Chart_Type.$TSV_File;
				$inFile = File_get_contents( $File );
				$inFile_Lines = Explode( "\n", $inFile );
				$Anno_Array = Array();
				$miRNA_Array = Array();

				For( $j = 1; $j < Count( $inFile_Lines )-1; ++$j )
				{
					$inFile_Line = Explode( "\t", $inFile_Lines[$j] );
					$anno = Explode( ':', $inFile_Line[0] );

					if( $isAbundant == 'MostAbundant' )
					{
						$miRNA = Explode( '*', $anno[0] );

						if( Count( $miRNA ) != 2 )
							Continue;
					}

					Array_Push( $Anno_Array, $anno[0] );

					if( $miRNA_Select == $anno[0] )
					{
						$Tail_Array = Array();

						for( $k = 1; $k < Count( $inFile_Line ); ++$k )
						{
							Array_Push( $Tail_Array, $inFile_Line[$k] );
						}

						$miRNA_Array[ $anno[1] ] = $Tail_Array;
					}
				}

				Sort( $Anno_Array, SORT_STRING );
				$Uniq_Anno_Array = Array_Unique( $Anno_Array, SORT_STRING );

				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';

				echo '<select name=miRNA_Select onchange=this.form.submit();>';
				echo '<option '; if($miRNA_Select=='') echo 'selected'; echo '>Select miRNA</option>';

				For( $i = 0; $i < Count( $Anno_Array ); ++$i )
				{
					if( $Uniq_Anno_Array[$i] != '' )
					{
						echo '<option value='.$Uniq_Anno_Array[$i].' ';

						if( $miRNA_Select == $Uniq_Anno_Array[$i] ) 
							echo 'selected ';

						echo '>'.$Uniq_Anno_Array[$i].'</option>';
					}
				}

				echo "</select>
					<input type='hidden' name='ForceY' value='$ForceY' />
					<input type='hidden' name='TSV_File' value='$TSV_File' />
					<input type='hidden' name='Chart_Type' value='$Chart_Type' />
					</form>";

#<!--================== ForceY ====================-->

				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';
				echo '<input type=text name=ForceY size=3 value=';

				if( $ForceY=='' )
					echo 'Hight';
				else
					echo $ForceY;

				echo " onfocus=\"{this.value='';}\">";

				echo "</select>
					<input type='hidden' name='TSV_File' value='$TSV_File' />
					<input type='hidden' name='Chart_Type' value='$Chart_Type' />
					<input type='hidden' name='miRNA_Select' value='$miRNA_Select' />
					<input type='submit' value='Submit' /> 
					</form>";

#<!--================== miRNA Tail Bar Chart ====================-->

				$index = Explode( "\t", $inFile_Lines[0] );
				$miRNA_Keys = Array_keys( $miRNA_Array );

				$Temp = Tempnam( '/tmp', $miRNA_Select );
				$Ftemp = fopen( $Temp, 'w' );
				Fwrite( $Ftemp, $index[0] );

				For( $i = 0; $i < Count( $miRNA_Keys ); $i++ )
				{
					Fwrite( $Ftemp, "\t".$miRNA_Keys[$i] );
				}

				Fwrite( $Ftemp, "\n".$index[6] );

				For( $i = 0; $i < Count( $miRNA_Keys ); $i++ )
				{
					Fwrite( $Ftemp, "\t".$miRNA_Array[ $miRNA_Keys[$i] ][5] );
				}

				Fwrite( $Ftemp, "\n" );

				For( $j = 1; $j < Count( $index )-1; $j++ )
				{
					Fwrite( $Ftemp, $index[$j] );

					For( $i = 0; $i < Count( $miRNA_Keys ); $i++ )
					{
						Fwrite( $Ftemp, "\t".$miRNA_Array[ $miRNA_Keys[$i] ][$j-1] );
					}

					Fwrite( $Ftemp, "\n" );
				}

				fclose( $Ftemp );

				echo '<svg id=bar></svg>';
				echo "<script>
						d3.tsv( '$Temp', function( tsv_data ) {
						var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'miRNA_Length'; });
						var data = new Array();
			
						tsv_data.forEach( function(d) {
							d.len = length.map( function(key) {
								return {
									x: Number(key),
									y: +d[key]
								};
			
							});
							var sample = new Object();
							sample.key = d['miRNA_Length'];
							sample.values = d.len;
							data.push( sample );
						});
			
						nv.addGraph({
							generate: function() {
								var width = nv.utils.windowSize().width,
									height = nv.utils.windowSize().height;
					
								var chart = nv.models.multiBarChart()";

					if( $ForceY != '' && $ForceY != 'Hight' )
						echo "		.forceY([$ForceY,0])";

				echo "				.width(width)
									.height(height)
									.color( ['#000000', '#FF0000', '#0000FF', '#FFBF00', '#088A08', '#6E6E6E'])
									.stacked(true);
				
								var svg = d3.select('#bar').datum(data);
								console.log('calling chart');
								svg.transition().duration(0).call(chart);
					
								return chart;
							},
							callback: function(graph) {
									nv.utils.windowResize(function() {
										var width = nv.utils.windowSize().width;
										var height = nv.utils.windowSize().height;
										graph.width(width).height(height);

										d3.select('#bar')
												.attr('width', width)
												.attr('height', height)
												.transition().duration(0)
												.call(graph);
									});
							}
						});
					});
				</script>";
			}
		}
	?>
	</body>
</html>
