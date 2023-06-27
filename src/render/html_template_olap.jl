const html_template_olap=mt"""<!DOCTYPE html>
<html>
<meta charset="utf-8">
<head>
    <title></title>
    <style>
    abbr[title] {
    text-decoration: none;
    }
    table, td {
        border-collapse: collapse;
        margin: 15px 15px;
        padding: 5px 5px;
        table-layout: fixed;
        min-width: 85px;
    }
    .top_row {
    font-weight: bold;
    color: #808080;
    }

    thead,tfoot {
        font-weight: bold;
        background-color: #333;
        color:white;
    }

    .info {
        background-color: #E2E2E2;
        margin:5px;
        padding:5px;
    }
    </style>
    <script id="MathJax-script" async
        src="https://cdn.jsdelivr.net/npm/mathjax@3.0.0/es5/tex-mml-chtml.js">
</script>
</head>

<body>
    <div style="display:flex;">
        <div style="float:left; margin:25px; border:1px solid black; max-width:3500px; padding:10px;" > 			
            Number of sequences: {{:num_seq}} <br>
            <table>
                <thead>
                    <tr>
                        <th colspan="100%">
                            Discovered motifs
                        </th>
                    </tr>
                </thead>
            <tbody>	
                    <tr class="top_row">
                        <td style="text-align: center"><abbr title="A label assigned to each discovered motif">Label</abbr></td>
                        <td style="text-align: center"><abbr title="Significance of the enrichment of each motif. Done via fisher exact test.">
                                        P-value</abbr></td>
                        <td style="text-align: center"><abbr title="An estimate of the number of instances in the dataset"># instances</abbr></td>
                        <td style="text-align: center"><abbr title="Position weight matrix">Logo</abbr></td>            
                    </tr>		

                    {{#:DF}}
                    <tr>
                        <td style="text-align:center"><a href="{{{:logo_folder}}}/{{:logo}}.transfac">{{:label}}</a></td>
                        <td>{{{:eval}}}</td>
                        <td style="text-align:center">{{{:counts}}}</td>
                        <td><img id="d_logo_{{:label}}" height=65 src="{{{:logo_folder}}}/{{:logo}}.png"><br>
                            <div id="d_orientation_{{:label}}"></div><br>
                            <button type="button" onclick="discovered_{{:label}}_changeToRC()">Reverse complement</button>
                        </td>                
                        <script type="text/javascript">					
                                    function discovered_{{:label}}_changeToRC() {
                                        var image = document.getElementById("d_logo_{{:label}}");
                                        if (image.src.match("_c")) {
                                            image.src = "{{{:logo_folder}}}/{{:logo}}.png";
                                        } else {
                                            image.src = "{{{:logo_folder}}}/{{:logo}}_c.png";
                                        }
                                        var orientation = document.getElementById("d_orientation_{{:label}}");
                                        if (orientation.innerHTML === ""){
                                            orientation.innerHTML = "reverse-complement";
                                        } else {
                                            orientation.innerHTML = "";
                                        }
                                    }
                        </script>
                    </tr>
                    {{/:DF}}
                </tbody>
            </table>
            <br><br>
        </div>
    </div>
</body>
</html>
"""