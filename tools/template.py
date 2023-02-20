import pandas as pd
from io import StringIO
from jinja2 import Template

def format_output(data, columns, lat):

    formated_data = None

    """
    format_output() function converts the output of any language or tool into the same format: Pandas dataframe.
    :param data: output data
    :param columns: column names
    :param lat: language or tool
    :return: a Pandas dataframe
    """

    if lat == "python" or lat == "cypher":
        formated_data = data
    elif lat == "sql":
        formated_data = pd.DataFrame(data=data, columns=columns)
    elif lat == "sparql":
        formated_data = pd.DataFrame(data).applymap(lambda x: x["value"].replace("http://MGDB.com/p", ""))
    elif lat == "clingo":
        formated_data = pd.read_csv(StringIO("\n".join([x.strip()[1:-1] for x in data])), names=columns)
    else:
        raise ValueError("lat parameter only supports [python, sql, cypher, sparql, clingo] for now.")

    return formated_data

def report_output_per_query(formated_data, title, edge, include_interactive_table, include_network_graph):

    base_html = """
    <!doctype html>
    <html lang="en">
      <head>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1">
        <title>{title}</title>
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css">
        <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/js/bootstrap.min.js"></script>
      </head>
      <body>
        <div class="container">
            {interactive_table_content}
            {network_graph_content}
        </div>
        {interactive_table_script}
        {network_graph_script}
      </body>
    </html>
    """

    (interactive_table_content, interactive_table_script) = generate_interactive_table(formated_data, title) if include_interactive_table else ("", "")
    (network_graph_content, network_graph_script) = generate_network_graph(formated_data, title, edge) if include_network_graph else ("", "")

    html = base_html.format(title=title, 
                            interactive_table_content=interactive_table_content,
                            network_graph_content=network_graph_content,
                            interactive_table_script=interactive_table_script,
                            network_graph_script=network_graph_script)
        
    return html

    
    
def generate_interactive_table(formated_data, title):

    """
    generate_interactive_table() function creates an interactive HTML table with a Pandas dataframe.
    :param formated_data: a Pandas dataframe
    :param title: title of the table
    :return interactive_table_content: the HTML code within the <div class="container"></div> tag
    :return interactive_table_script: javascripts (e.g., bootstrap-table.js) and css files used to generate an interactive table
    """

    base_interactive_table_content = """
            <h3 style="color: #2E6F9F">An interactive table of {{title}}</h3>
            <div id="toolbar">
                <select class="form-control">
                        <option value="">Export Data</option>
                        <option value="all">Export All</option>
                        <option value="selected">Export Selected</option>
                </select>
            </div>
            <table id="table" data-toggle="table" data-height="460" data-search="true" data-filter-control="true" data-show-export="true" data-click-to-select="true" data-toolbar="#toolbar" class="table-responsive">
                <thead>
                    <tr>
                        <th data-checkbox="true"></th>
                        {% for col_name in df.columns %}
                            <th data-field="{{col_name}}" data-filter-control="select" data-sortable="true">{{ col_name }}</th>
                        {% endfor %}
                    </tr>
                </thead>
                <tbody>
                    {% for idx, row in df.iterrows() %}
                    <tr>
                        <td class="bs-checkbox "><input data-index="{{ idx }}" name="btSelectItem" type="checkbox"></td>
                        {% for col_name in df.columns %}
                        <td>{{ row[col_name] }}</td>
                        {% endfor %}
                    </tr>
                    {% endfor %}
                </tbody>
            </table>
    """

    interactive_table_script = """
        <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.10.0/bootstrap-table.min.css">
        <link rel="stylesheet" href="https://rawgit.com/vitalets/x-editable/master/dist/bootstrap3-editable/css/bootstrap-editable.css">
        <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/2.1.3/jquery.min.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.10.0/bootstrap-table.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.9.1/extensions/editable/bootstrap-table-editable.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.9.1/extensions/export/bootstrap-table-export.js"></script>
        <script src="http://frontendfreecode.com/codes/files/tableExport.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.9.1/extensions/filter-control/bootstrap-table-filter-control.js"></script>
        <script id="rendered-js">
            var $table = $("#table");
            $(function () {
                $("#toolbar")
                    .find("select")
                    .change(function () {
                        $table.bootstrapTable("refreshOptions", {
                            exportDataType: $(this).val()
                        });
                    });
                });

                var trBoldBlue = $("table");

                $(trBoldBlue).on("click", "tr", function () {
                    $(this).toggleClass("bold-blue");
                });
        </script>
    """

    tpl = Template(base_interactive_table_content)
    interactive_table_content = tpl.render(df=formated_data, title=title)

    return (interactive_table_content, interactive_table_script)

def generate_network_graph(formated_data, title, edge):

    """
    generate_network_graph() function visualize a network graph with a Pandas dataframe.
    :param formated_data: a Pandas dataframe
    :param title: title of the table
    :param edges: a list of list. [[start_node1, end_node1], [start_node2, end_node2], ...]
    :return network_graph_content: the HTML code within the <div class="container"></div> tag
    :return network_graph_script: javascripts (e.g., highcharts.js) used to visualize the network graph
    """

    if edge == []:
        return ('<br><h3 style="color: #2E6F9F">No network graph was generated due to the lack of edges.</h3>', '')

    network_graph_content = """
    <br><br><h3 style="color: #2E6F9F">A Force-Directed Network Graph of {title}</h3>
    <figure class="highcharts-figure">
        <div id="graph-container"></div>
        <p class="highcharts-description"></p>
    </figure>
    """.format(title=title)

    base_network_graph_script = """
    <script src="https://code.highcharts.com/highcharts.js"></script>
    <script src="https://code.highcharts.com/modules/networkgraph.js"></script>
    <script src="https://code.highcharts.com/modules/exporting.js"></script>
    <script src="https://code.highcharts.com/modules/accessibility.js"></script>
    <script>
        Highcharts.addEvent(
            Highcharts.Series,
            'afterSetOptions',
            function (e) {
                var colors = Highcharts.getOptions().colors,
                    i = 0,
                    nodes = {},
                    l2r = {};
                    r2l = {};
                if (
                    this instanceof Highcharts.Series.types.networkgraph &&
                    e.options.id === 'network_graph'
                ) {
                    e.options.data.forEach(function (link) {
                        if (link[0] in l2r) {
                            nodes[link[1]] = {
                                id: link[1],
                                color: nodes[l2r[link[0]]].color
                            };
                        } else if (link[1] in r2l) {
                            nodes[link[0]] = {
                                id: link[0],
                                color: nodes[r2l[link[1]]].color
                            };
                        } else {
                            l2r[link[0]] = link[1];
                            r2l[link[1]] = link[0];
                            nodes[link[0]] = {
                                id: link[0],
                                color: colors[i]
                            };
                            i = (i + 1) % 10;
                            nodes[link[1]] = {
                                id: link[1],
                                color: colors[i++]
                            };
                        }

                    });

                    e.options.nodes = Object.keys(nodes).map(function (id) {
                        return nodes[id];
                    });
                }
            }
        );

        Highcharts.chart('graph-container', {
            chart: {
                type: 'networkgraph',
                height: '100%'
            },
            title: {text: undefined},
            plotOptions: {
                networkgraph: {
                    keys: ['from', 'to'],
                    layoutAlgorithm: {
                        enableSimulation: true,
                        friction: -0.9
                    }
                }
            },
            series: [{
                accessibility: {
                    enabled: false
                },
                dataLabels: {
                    enabled: true,
                    linkFormat: ''
                },
                id: 'network_graph',
                data: [
                    {% for e in edge %}
                        {% for _, row in df.iterrows() %}
                        ["{{row[e[0]]}}","{{row[e[1]]}}"],
                        {% endfor %}
                    {% endfor %}
                ]
            }]
        });
    </script>
    """
    tpl = Template(base_network_graph_script)
    network_graph_script = tpl.render(df=formated_data, edge=edge)

    return (network_graph_content, network_graph_script)
