
//// This document contains the JavaScript code for the web application. It is quite large, and should be split into multiple sections for better readability and maintainability. ////

// Fetch data from data.json - this should be done before the rest of the script runs
fetch("/static/data.json")
    .then(response => response.json())
    .then(data => {
        window.query_targets = data.query_targets;
        window.qt_pairs = data.qt_pairs;
        window.nodes = data.nodes;
        window.node_clusters = data.node_clusters;
        window.afdb_hits = data.afdb_hits;
        window.vzone_hits = data.vzone_hits;
        window.species_colors = data.species_colors;
        window.query_superposed_pdb = data.query_superposed_pdb;
        window.vzone_heatmap = data.vzone_heatmap;
        window.cono_heatmap = data.cono_heatmap;
        window.cono_hits = data.cono_hits;
    })
    .catch(err => console.error("Error loading data.json:", err));

    /// SIDE-PANEL RESIZER ///
    // Resize functionality for the side panel - allows users to adjust the width of the side panel
    const resizer = document.getElementById("resizer");
    const sidePanel = document.getElementById("side-panel");
    const mainContainer = document.getElementById("main-container");

    let isResizing = false;
    let startX = 0;
    let startWidth = 0;

    resizer.addEventListener("mousedown", function (e) { // Start resizing on mousedown
    isResizing = true;
    startX = e.clientX;
    startWidth = parseInt(window.getComputedStyle(sidePanel).width, 10);
    document.body.style.cursor = "col-resize";
    e.preventDefault(); // Prevent text selection while resizing
    });

    document.addEventListener("mousemove", function (e) {
    if (!isResizing) return;

    const dx = startX - e.clientX; // Calculate the change in mouse position
    const newWidth = startWidth + dx; // Calculate the new width of the side panel

    if (newWidth > 250 && newWidth < window.innerWidth * 0.7) {
        sidePanel.style.width = newWidth + "px";
    }
    });

    document.addEventListener("mouseup", function () { // Stop resizing on mouseup
    if (isResizing) {
        isResizing = false;
        document.body.style.cursor = "default";
        }
    });


        /// RENDER PDBs ///
        // This viewer will be used to display PDB structures
        var currentNode = null;

        function loadPDBFromFile(nodeName) {
        let url = `/get_pdb/${encodeURIComponent(nodeName)}`;
        return fetch(url) // Fetch the PDB file for the given node name
            .then(response => {
            if (!response.ok) {
                throw new Error("PDB file not found for node: " + nodeName);
            }
            return response.text();
            });
        }

        function loadSuperposedPDB(key) {
        let url = window.query_superposed_pdb[key];
        if (!url) {
            return Promise.reject("No superposed PDB for key: " + key);
        }
        return fetch(url) // Fetch the superposed PDB file
            .then(response => {
            if (!response.ok) throw new Error("File not found: " + url);
            return response.text();
            });
        }

        // Render the PDB structures in the 3Dmol.js viewer 
        function renderPDB(finalQueryCoords, finalTargetCoords, nodeName, targetList) {
            let viewerContainer = document.getElementById('structure-viewer');
            viewerContainer.innerHTML = "";
            let viewer = $3Dmol.createViewer("structure-viewer", {
                defaultcolors: $3Dmol.rasmolElementColors
            });
            document.getElementById('structure-viewer').style.top = 'auto';
            document.getElementById('structure-viewer').style.left = 'auto';
            viewer.removeAllModels();
            try { // Attempt to add the models to the viewer
                let queryModel = viewer.addModel(finalQueryCoords, "pdb");
                queryModel.setStyle({cartoon: {color: '#d80032'}});
                let targetModel = viewer.addModel(finalTargetCoords, "pdb");
                targetModel.setStyle({cartoon: {color: '#ee9b00'}});
            } catch (error) {
                console.error("Error loading models in 3Dmol.js", error);
            }
            viewer.zoomTo(); // Zoom to fit the models in the viewer
            viewer.render(); // Render the viewer

            let key = "";
            if (targetList.length > 0) {
                key = nodeName + "__" + targetList[0];
            } else {
                key = nodeName + "||" + nodeName;
            }

            updatePhyloConnections(nodeName); // Update phylogenetic connections based on the current node
        }

        function extractUniprotId(target) { // Extract the UniProt ID from the target string
            var parts = target.split("-");
            return (parts[0] === "AF" && parts.length >= 2) ? parts[1] : target;
        }

        function updateSidePanel(nodeName) { // Update the side panel with information about the selected node
            currentNode = nodeName;
            document.getElementById('node-name').innerHTML = "Node: " + nodeName;

            // Render content in the side panel
            renderNodeClassification(nodeName);
            renderClusterClassification(nodeName);
            loadAFDBDetailsForNode(nodeName);
            renderConotoxinNodeClassification(nodeName);
            renderConotoxinClusterClassification(nodeName);

            // Update the structure viewer with the current node and its targets
            let targetList = window.query_targets[nodeName] || [];

            if (!window.currentTargetIndex) {
                window.currentTargetIndex = {};
            }
            if (window.currentTargetIndex[nodeName] === undefined) {
                window.currentTargetIndex[nodeName] = 0;
            }
            let currentIndex = window.currentTargetIndex[nodeName];

            // Update the target index display - shows the current target and total number of targets
            document.getElementById("target-index").innerHTML = targetList.length > 0 ? (currentIndex + 1) + " / " + targetList.length : "";

            document.getElementById("node-source").innerHTML = "Node (<span style=color:#d80032;>&#9679;</span>): " + shortName(nodeName);
            if (targetList.length > 0) {
                document.getElementById("target-source").innerHTML = "Target (<span style=color:#ee9b00;>&#9679;</span>): " + shortName(targetList[currentIndex]);
            } else {
                document.getElementById("target-source").innerHTML = "Target (<span style=color:#ee9b00;>&#9679;</span>): " + shortName(nodeName);
            }

            // Load the PDB files for the current node and its target
            loadPDBFromFile(nodeName)
            .then(queryPDB => {
                if (targetList.length > 0) {
                    let targetId = targetList[currentIndex];
                    return Promise.all([queryPDB, loadPDBFromFile(targetId)]);
                }
                return Promise.resolve([queryPDB, queryPDB]);
            })
            .then(([queryPDB, targetPDB]) => {
                let key = (targetList.length > 0) ? nodeName + "||" + targetList[currentIndex] : nodeName + "||" + nodeName;
                loadSuperposedPDB(key)
                .then(superposedPDB => {
                    renderPDB(superposedPDB, targetPDB, nodeName, targetList);
                })
                .catch(error => {
                        console.error("Error loading superposed PDB:", error);
                        renderPDB(queryPDB, targetPDB, nodeName, targetList);
                });
            })
            .catch(error => {
                console.error("Error loading PDB data:", error);
                document.getElementById('structure-viewer').innerHTML = "<p>Error: Missing structure data.</p>";
            });
        }

        // Switch targets in Structure Viewer
        // These functions update the structure viewer with the next or previous target
        document.getElementById("prev-target").onclick = function() {
            if (currentNode) {
                let targetList = window.query_targets[currentNode] || [];
                if (targetList.length > 0) {
                    if (!window.currentTargetIndex) {
                        window.currentTargetIndex = {};
                    }
                    if (window.currentTargetIndex[currentNode] === undefined) {
                        window.currentTargetIndex[currentNode] = 0;
                    }
                    let index = window.currentTargetIndex[currentNode];
                    index = (index - 1 + targetList.length) % targetList.length;
                    window.currentTargetIndex[currentNode] = index;
                    updateSidePanel(currentNode);
                }
            }
        };

        document.getElementById("next-target").onclick = function() {
            if (currentNode) {
                let targetList = window.query_targets[currentNode] || [];
                if (targetList.length > 0) {
                    if (!window.currentTargetIndex) {
                        window.currentTargetIndex = {};
                    }
                    if (window.currentTargetIndex[currentNode] === undefined) {
                        window.currentTargetIndex[currentNode] = 0;
                    }
                    let index = window.currentTargetIndex[currentNode];
                    index = (index + 1) % targetList.length;
                    window.currentTargetIndex[currentNode] = index;
                    updateSidePanel(currentNode);
                }
            }
        };

        // Helper function to extract a short name from a longer node name
        function shortName(name) {
            var match = name.match(/(test_[a-zA-Z]+_\\d+)/);
            return match ? match[1] : name;
        }

        // Function to get the color for the AF heatmap based on the E-value - subjectively selected thresholds
        function getColorAFHeatmap(e) {
            if (e == null) {
                return "#ffffff";
            }
            const logVal = -Math.log10(e);
            if (logVal >= 50) {
                return "#0d47a1";
            }
            if (logVal >= 40) {
                return "#1976d2";
            }
            else if (logVal >= 30) {
                return "#2196f3";
            }
            else if (logVal >= 20) {
                return "#64b5f6";
            }
            else if (logVal >= 10) {
                return "#90caf9";
            }
            else if (logVal >= 5) {
                return "#bbdefb";
            }
            else {
                return "#e3f2fd";
            }
        }

        function getColorConoHeatmap(e) {
            if (e == null) {
                return "#ffffff";
            }
            const logVal = -Math.log10(e);
            if (logVal >= 50) {
                return "#155d27";
            }
            else if (logVal >= 40) {
                return "#208b3a";
            }
            else if (logVal >= 30) {
                return "#25a244";
            }
            else if (logVal >= 20) {
                return "#4ad66d";
            }
            else if (logVal >= 10) {
                return "#6ede8a";
            }
            else if (logVal >= 5) {
                return "#92e6a7";
            }
            else {
                return "#b7efc5";
            }
        }


        // Render the AF heatmap using D3.js
        function renderHeatmap(heatmapData) {
                // Ensure data is present
                heatmapData = heatmapData || window.vzone_heatmap;
                console.log("D3 heatmap data:", heatmapData); // Debug logging
                if (!heatmapData || Object.keys(heatmapData).length === 0) {
                document.getElementById("af-heatmap-content").innerHTML = "<p>No heatmap data available.</p>"; // Text displayed if no data is available
                return;
                }
            
                const species = Object.keys(window.species_colors);
                const classifications = Object.keys(heatmapData);
                
                // Prepare a flat array for D3
                const data = [];
                species.forEach(function(sp) {
                    classifications.forEach(function(cl) {
                    data.push({
                        species: sp,
                        classification: cl,
                        value: heatmapData[cl][sp]
                    });
                });
            });
        
            // Define dimensions and margins
            const margin = { top: 200, right: 100, bottom: 120, left: 150 },
                    cellSize = 30,
                    width = classifications.length * cellSize,
                    height = species.length * cellSize;
            
            // Remove any existing SVG from the container
            d3.select("#af-heatmap-content").selectAll("svg").remove();
            
            // Append an SVG element
            const svg = d3.select("#af-heatmap-content")
                            .append("svg")
                            .attr("width", width + margin.left + margin.right)
                            .attr("height", height + margin.top + margin.bottom)
                            .append("g")
                            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
            
            // Create scales for the x and y axes
            const x = d3.scaleBand()
                        .domain(classifications)
                        .range([0, width])
                        .padding(0.05);
            
            const y = d3.scaleBand()
                        .domain(species)
                        .range([0, height])
                        .padding(0.05);
            
            // Draw grid cells
            svg.selectAll("rect")
                .data(data)
                .enter().append("rect")
                .attr("x", d => x(d.classification))
                .attr("y", d => y(d.species))
                .attr("width", x.bandwidth())
                .attr("height", y.bandwidth())
                .style("fill", d => getColorAFHeatmap(d.value))
                .style("stroke", "#fff");
            
            const tooltip = d3.select("body")
                            .append("div")
                                .attr("class", "heatmap-tooltip")
                                .style("position", "absolute")
                                .style("z-index", "10")
                                .style("visibility", "hidden")
                                .style("background", "#f8f8f8")
                                .style("padding", "5px")
                                .style("border", "1px solid #ccc")
                                .style("border-radius", "3px");
            
            svg.selectAll("rect")
                .on("mouseover", function(event, d) {
                    tooltip.style("visibility", "visible")
                        .text(`Species: ${d.species}, Classification: ${d.classification}, E-value: ${d.value != null ? d.value : "N/A"}`);
                })
                .on("mousemove", function(event) {
                    tooltip.style("top", (event.pageY - 10) + "px")
                        .style("left", (event.pageX + 10) + "px");
                })
                .on("mouseout", function() { tooltip.style("visibility", "hidden"); });
            
            svg.append("g")
                .attr("class", "x axis")
                .call(d3.axisTop(x))
                .selectAll("text")
                .style("text-anchor", "start")
                .attr("transform", "rotate(-45)")
                .attr("dx", "0.5em")
                .attr("dy", "-0.5em");
            
            // Species specific images to display along y-axis
            // Need to be present in images/ directory, and added along with the identifier
            const speciesImages = {
            "test_crotalus": "images/crotalus.png",
            "test_sponsalis": "images/sponsalis.png",
            "test_sting": "images/neotrygon.png",
            "test_vespa": "images/vespa.png",
            "test_shrew": "images/blarina.png",
            "test_latrohesperus": "images/latrodectus.png",
            "test_najanaja": "images/najanaja.png",
            "test_heloderma": "images/heloderma.png",
            "test_myrmica": "images/myrmica.png",
            "test_geographus": "images/geographus.png",
            "test_blenny": "images/meiacanthus.png",
            "test_geographus": "images/geographus.png",
            "test_scolopendra": "images/scolopendra.png",
            "test_platymeris": "images/platymeris.png",
            "test_euscorpius": "images/euscorpius.png",
            "test_gloriamaris": "images/gloriamaris.png",
            "test_xibalbanus": "images/xibalbanus.png",
            "test_unedogemmula": "images/unedogemmula.png",
            "test_gemmula": "images/gemmula.png",
            "test_hapalochlaena": "images/hapalochlaena.png",
            "test_crassispira": "images/crassispira.png",
            "test_actinia": "images/actinia.png",
            "test_corixa": "images/corixa.png"
            };

            // Add y-axis (speciesImages) on left
            const yAxis = svg.append("g")
                            .attr("class", "y axis")
                            .call(d3.axisLeft(y));

            yAxis.selectAll(".tick")
                .each(function(speciesName) {
                const imgUrl = speciesImages[speciesName] || "https://via.placeholder.com/30";
                d3.select(this)
                    .append("svg:image")
                    .attr("xlink:href", imgUrl)
                    .attr("x", -130)  // adjust offset to position image left of the text
                    .attr("y", - (y.bandwidth()) / 2)  // center vertically within the band
                    .attr("width", 30)
                    .attr("height", 30);
                });
        }

        // Render the conotoxin heatmap using D3.js
        function renderConoHeatmapD3(heatmapData) {
            heatmapData = heatmapData || window.cono_heatmap;
            if (!heatmapData || Object.keys(heatmapData).length === 0) {
                document.getElementById("cono-heatmap-content")
                        .innerHTML = "<p>No conotoxin data available.</p>";
                return;
            }

            // rows = species, cols = conotoxin superfamilies
            const species         = Object.keys(window.species_colors);
            const superfamilies   = Object.keys(heatmapData);
            const data = [];
            species.forEach(sp => {
                superfamilies.forEach(sf => {
                data.push({
                    species:       sp,
                    classification: sf,
                    value:         heatmapData[sf][sp]
                });
                });
            });

            const margin = { top: 140, right: 100, bottom: 120, left: 150 },
                    cellSize = 30,
                    width    = superfamilies.length * cellSize,
                    height   = species.length       * cellSize;

            d3.select("#cono-heatmap-content").selectAll("svg").remove();

            const svg = d3.select("#cono-heatmap-content")
                .append("svg")
                .attr("width",  width + margin.left + margin.right)
                .attr("height", height + margin.top  + margin.bottom)
                .append("g")
                .attr("transform", `translate(${margin.left},${margin.top})`);

            const x = d3.scaleBand()
                        .domain(superfamilies)
                        .range([0, width])
                        .padding(0.05);
            const y = d3.scaleBand()
                        .domain(species)
                        .range([0, height])
                        .padding(0.05);

            svg.selectAll("rect")
                .data(data)
                .enter().append("rect")
                .attr("x", d => x(d.classification))
                .attr("y", d => y(d.species))
                .attr("width",  x.bandwidth())
                .attr("height", y.bandwidth())
                .style("fill", d => getColorConoHeatmap(d.value))
                .style("stroke", "#FFF")
                .on("mouseover", function(event, d){
            });

            svg.append("g")
                .call(d3.axisTop(x))
                .selectAll("text")
                .attr("transform", "rotate(-45)")
                .style("text-anchor", "start");
            svg.append("g")
                .call(d3.axisLeft(y));
            }

            window.onload = function() {
                let network = window.network;
                network.on("click", function(properties) {
                    if (properties.nodes.length > 0) {
                        let clickedNode = properties.nodes[0];
                        if(!window.currentTargetIndex) window.currentTargetIndex = {};
                        window.currentTargetIndex[clickedNode] = 0;
                        updateSidePanel(clickedNode);
                    }
                });
        };

        window.network.once("stabilizationIterationsDone", function() {
            window.network.setOptions({ physics: false });
        });

        let connectionGroup = null;
        
        function parseNewick(newick) {
        // Remove any trailing semicolon and trim whitespace
        newick = newick.trim();
        if (newick.endsWith(";")) {
            newick = newick.slice(0, -1);
        }
        // Split the string into tokens
        const tokens = newick.split(/\s*(;|\(|\)|,|:)\s*/).filter(t => t && t !== "");
        let ancestors = [];
        let tree = {};
        
        // Process each token
        for (let i = 0; i < tokens.length; i++) {
            const token = tokens[i];
            switch (token) {
            case "(":
                // Start a new subtree
                let subtree = {};
                if (!tree.children) {
                tree.children = [];
                }
                tree.children.push(subtree);
                ancestors.push(tree);
                tree = subtree;
                break;
            case ",":
                // Start a new sibling
                let sibling = {};
                ancestors[ancestors.length - 1].children.push(sibling);
                tree = sibling;
                break;
            case ")":
                // End the current subtree
                tree = ancestors.pop();
                break;
            case ":":
                // Skip branch lengths
                i++; // jump over the branch length token.
                break;
            default:
                tree.name = token.replace(/^'+/, "").replace(/'+$/, "");
            }
        }
        return tree;
        }

        // Species latin names and their corresponding identifiers
        const speciesMapping = {
            'Conus sponsalis': 'test_sponsalis',
            'Crotalus oreganus': 'test_crotalus',
            'Neotrygon kuhlii': 'test_sting',
            'Vespa velutina': 'test_vespa',
            'Blarina brevicauda': 'test_shrew',
            'Naja naja': 'test_najanaja',
            'Latrodectus hesperus': 'test_latrohesperus',
            'Heloderma horridum': 'test_heloderma',
            'Myrmica rubra': 'test_myrmica',
            'Conus geographus': 'test_geographus',
            'Meiacanthus atrodorsalis': 'test_blenny',
            'Euscorpius italicus': 'test_euscorpius',
            'Scolopendra subspinipes': 'test_scolopendra',
            'Platymeris biguttatus': 'test_platymeris',
            'Conus gloriamaris': 'test_gloriamaris',
            'Corixa punctata': 'test_corixa',
            'Unedogemmula bisaya': 'test_unedogemmula',
            'Gemmula speciosa': 'test_gemmula',
            'Crassispira cerithina': 'test_crassispira',
            'Xibalbanus tulumensis': 'test_xibalbanus',
            'Hapalochlaena maculosa': 'test_hapalochlaena',
            'Actinia tenebrosa': 'test_actinia'
        };

        // Parse the Newick string
        const newickString = "(((('Platymeris biguttatus':4,'Corixa punctata':4)Hemiptera:4,('Vespa velutina':4,'Myrmica rubra':4)Hymenoptera:4)Insecta:4,('Latrodectus hesperus':4,'Euscorpius italicus':4)Arachnida:4,(('Scolopendra subspinipes':4)Chilopoda:4),('Xibalbanus tulumensis':4)Remipedia:4)Arthropoda:4,(('Meiacanthus atrodorsalis':4)Actinopteri:4,('Neotrygon kuhlii':4)Chondrichthyes:4,('Blarina brevicauda':4)Mammalia:4,(('Crotalus oreganus':4)Viperidae:4,('Naja naja':4)Elapidae:4,('Heloderma horridum':4)Helodermatidae:4)Squamata:4)Chordata:4,(('Hapalochlaena maculosa')Cephalopoda:4,(('Crassispira cerithina':4,'Unedogemmula bisaya':4,'Gemmula speciosa':4)Turridae:4,('Conus gloriamaris':4,'Conus sponsalis':4,'Conus geographus':4)Conus:4)Neogastropoda:4)Mollusca:4,(('Actinia tenebrosa':4)Cnidaria:4)Metazoa:4);";
        const treeData = parseNewick(newickString);

        // Create a d3 cluster layout for a dendrogram
        const margin = {top: 20, right: 20, bottom: 20, left: 20};

        // Define a desired ultrametric length for branch lengths
        const ultrametricLength = 500;  

        // Set the SVG width to accommodate that length
        const numLeaves = d3.hierarchy(treeData).leaves().length;
        const rowHeight = 30;
        const svgWidth = 800 + margin.left + margin.right;
        const svgHeight = numLeaves * rowHeight + margin.top + margin.bottom;

        const svg = d3.select("#phylo-tree").append("svg")
            .attr("width", svgWidth)
            .attr("height", svgHeight)
            .append("g")
            .attr("transform", `translate(${margin.left},${margin.top})`);

        // Use d3.cluster to create a dendrogram layout
        const clusterLayout = d3.cluster().size([svgHeight - margin.top - margin.bottom, ultrametricLength]);
        const root = d3.hierarchy(treeData);
        clusterLayout(root);

        // Compute maximum depth
        const maxDepth = d3.max(root.descendants(), d => d.depth);
        root.leaves().forEach(leaf => { leaf.depth = maxDepth; });
        root.descendants().forEach(d => {
            d.y = (d.depth / maxDepth) * ultrametricLength;
        });

        // Draw red connectors
        connectionGroup = svg.append("g").attr("class", "phylo-connections");
        root.descendants().forEach(d => {
            if (speciesMapping[d.data.name]) {
                d.data.testName = speciesMapping[d.data.name];
            }
        });

        // Draw links as straight lines
        const linkGenerator = d3.linkHorizontal()
            .x(d => d.y)
            .y(d => d.x);

        svg.selectAll(".link")
            .data(root.links())
            .enter()
            .append("path")
            .attr("class", "link")
            .attr("d", linkGenerator)
            .attr("fill", "none")
            .attr("stroke", "#ccc")
            .attr("stroke-width", 2)
            .attr("stroke-dasharray", "4,2");

        // Draw nodes as circles
        const nodeGroup = svg.selectAll(".node")
            .data(root.descendants())
            .enter()
            .append("g")
            .attr("class", "node")
            .attr("transform", d => `translate(${d.y},${d.x})`);

        nodeGroup.append("circle")
            .attr("r", 8)
            .attr("fill", d => {
                if (d.data.testName && speciesColors[d.data.testName]) {
                    return speciesColors[d.data.testName];
                }
                return "#fff";
            })
            .attr("stroke", "#007BFF")
            .attr("stroke-width", 2)
            .style("cursor", "pointer")
            .on("click", function(event, d) {
                const selectedSpecies = d.data.testName || d.data.name;
                updateNetworkColors(selectedSpecies);
            });

        // Add labels to nodes
        nodeGroup.append("text")
            .attr("dy", 4)
            .attr("x", d => d.children ? -17 : 17)
            .style("text-anchor", d => d.children ? "end" : "start")
            .style("font-size", "10px")
            .text(d => d.data.testName ? `${d.data.name} (${d.data.testName})` : d.data.name);

        // Update network colors function
            function updateNetworkColors(selectedSpecies) {
                const allNodes = window.network.body.data.nodes.get();
                const updates = allNodes.map(n => ({
                    id: n.id,
                    color: (n.species === selectedSpecies) ? (window.species_colors[n.species] || n.color) : 'gray'
                }));
                window.network.body.data.nodes.update(updates);
            }
        
        // Reset button to restore original colors
        document.getElementById("reset-species").onclick = function() {
            const allNodes = window.network.body.data.nodes.get();
            const updates = allNodes.map(n => ({
                id: n.id,
                color: speciesColors[n.species] || n.color
            }));
            window.network.body.data.nodes.update(updates);
            // Remove any drawn connection lines in the tree.
            connectionGroup.selectAll("line").remove();
        };

        root.each(d => {
        d.y = (d.depth / maxDepth) * (ultrametricLength*0.7);
        });

        let speciesPositions = {};
        root.descendants().forEach(d => {
        if (d.data.testName) {
            speciesPositions[d.data.testName] = { x: d.y, y: d.x };
            }
        });

        // Define a fixed x-value for the connector
        const maxX = d3.max(root.descendants(), d => d.y);
        const connectorX = maxX + 400;  // adjust the x value as needed

        function updatePhyloConnections(selectedNodeId) {
        // Remove previous connectors
        connectionGroup.selectAll("path").remove();
        
        // Get the species for the selected clustering node
        let selectedSpecies = window.network.body.data.nodes.get(selectedNodeId).species;
        
        // Gather the species of connected nodes
        let connectedIds = window.network.getConnectedNodes(selectedNodeId);
        let speciesSet = new Set();
        speciesSet.add(selectedSpecies);
        connectedIds.forEach(id => {
            let node = window.network.body.data.nodes.get(id);
            if (node && node.species) speciesSet.add(node.species);
        });
        let speciesArray = Array.from(speciesSet);
        
        // For each connected species, draw a connector from the selected species to that species
        speciesArray.forEach(species => {
            if (species === selectedSpecies) return;
            
            // Ensure both positions are available
            if (speciesPositions[selectedSpecies] && speciesPositions[species]) {
            let posA = speciesPositions[selectedSpecies];  // starting point (selected node species)
            let posB = speciesPositions[species];            // ending point (connected node species)
            
            // Build the path
            let pathData = `M ${posA.x + 370} ${posA.y} H ${connectorX} V ${posB.y} H ${posB.x + 370}`;
            
            connectionGroup.append("path")
                .attr("d", pathData)
                .attr("stroke", "red")
                .attr("stroke-width", 2)
                .attr("fill", "none");
            }
        });
        }

        // Navigation bar functionality
        document.querySelectorAll("#nav-bar li a").forEach(item => {
            item.addEventListener("click", function(e) {
                e.preventDefault();
                let target = this.getAttribute("href");
                // Hide all panels
                document.getElementById("node-details").style.display = "none";
                document.getElementById("phylo-tree-panel").style.display = "none";
                document.getElementById("af-heatmap-panel").style.display = "none";
                document.getElementById("classification-panel").style.display = "none";
                document.getElementById("upload-panel").style.display = "none";
                document.getElementById("export-panel").style.display = "none";

                // Show the selected panel based on the clicked link
                if(target === "#node"){
                    document.getElementById("node-details").style.display = "block";
                } else if(target === "#tree"){
                    document.getElementById("phylo-tree-panel").style.display = "block";
                } else if(target === "#af_classification"){
                    document.getElementById("classification-panel").style.display = "block";
                } else if(target === "#heatmap"){
                    document.getElementById("af-heatmap-panel").style.display = "block";
                    renderHeatmap(window.vzone_heatmap);
                    renderConoHeatmapD3(window.cono_heatmap);
                } else if(target === "#upload"){
                    document.getElementById("upload-panel").style.display = "block";
                } else if(target === "#export") {
                    document.getElementById("export-panel").style.display = "block";
                }
                document.querySelectorAll("#nav-bar li a").forEach(a => a.classList.remove("active"));
                this.classList.add("active");
            });
        });
        
        // Function to load AFDB details for a specific node
        function loadAFDBDetailsForNode(nodeName) {
            var rawEntries = window.afdb_hits && window.afdb_hits[nodeName] ? window.afdb_hits[nodeName] : [];
            rawEntries.sort(function(a, b) { return a.evalue - b.evalue; });
            var tableBody = document.querySelector("#afdb-table tbody");
            tableBody.innerHTML = "";  // Clear previous rows
            var limit = 5; // Limit to 5 entries - can be adjusted as needed
            rawEntries.slice(0, limit).forEach(function(entry) {
                var uniprotId = extractUniprotId(entry.target);
                fetch("https://rest.uniprot.org/uniprotkb/" + uniprotId + ".json") // Fetch UniProt data.json for entry
                .then(function(response) { return response.json(); })
                .then(function(data) { // Sometimes UniProt descriptors are not availabe - thus yielding N/A.
                    var name = (data.proteinDescription && data.proteinDescription.recommendedName &&
                            data.proteinDescription.recommendedName.fullName &&
                            data.proteinDescription.recommendedName.fullName.value) || uniprotId;
                    var mf = "N/A";
                    if (data.keywords) { // Check if Molecular GO-term is available
                        var mfKeywords = data.keywords.filter(function(kw) { return kw.category === "Molecular function"; });
                        if (mfKeywords.length > 0) {
                            mf = mfKeywords.map(function(kw) { return kw.name; }).join(", ");
                        }
                    }
                    // Append a new row to the table body with the fetched data
                    tableBody.innerHTML += "<tr><td>" + uniprotId + "</td><td>" + name + "</td><td>" + mf + "</td><td>" + entry.evalue + "</td></tr>";
                })
                .catch(function(error) {
                    console.error("Error fetching Uniprot data for " + uniprotId, error);
                    tableBody.innerHTML += "<tr><td>" + uniprotId + "</td><td>Error</td><td>Error</td><td>" + entry.evalue + "</td></tr>";
                });
            });
        }

        // Render VenomZone classification for a specific node
        function renderNodeClassification(nodeName) {
            var afdbMatches = window.afdb_hits && window.afdb_hits[nodeName] ? window.afdb_hits[nodeName] : [];
            afdbMatches.sort(function(a, b) { return a.evalue - b.evalue; });
            var afdbRows = "";
            afdbMatches.slice(0, 5).forEach(function(item) {
                afdbRows += "<tr><td>" + item.uniprot_id + "</td><td>" + item.name + "</td><td>" + item.molecular_function + "</td><td>" + item.evalue + "</td></tr>";
            });
            document.querySelector("#afdb-table tbody").innerHTML = afdbRows;
            
            var venomzoneMatches = window.vzone_hits && window.vzone_hits[nodeName] ? window.vzone_hits[nodeName] : [];
            venomzoneMatches.sort(function(a, b) { return a.evalue - b.evalue; });
            var venomzoneRows = "";
            venomzoneMatches.forEach(function(item, index) {
                if (index < 5) {
                    venomzoneRows += `<tr>
                                        <td>${item.classification || "N/A"}</td>
                                        <td>${item.evalue}</td>
                                    </tr>`;
                } else {
                    venomzoneRows += `<tr class="hidden-row" style="display:none;">
                                        <td>${item.classification || "N/A"}</td>
                                        <td>${item.evalue}</td>
                                    </tr>`;
                }
            });
            if (venomzoneMatches.length > 5) {
                venomzoneRows += `<tr class="show-more-row" id="venomzone-toggle" onclick="toggleVenomzoneRows()">
                                    <td colspan="2">Show more...</td>
                                </tr>`;
            }
            document.querySelector("#venomzone-table tbody").innerHTML = venomzoneRows;
        }

        // Toggle function for VenomZone rows - extend the table to show more classifications
        function toggleVenomzoneRows() {
            var extraRows = document.querySelectorAll("#venomzone-table tbody .hidden-row");
            var toggleRow = document.getElementById("venomzone-toggle");
            if (window.venomzoneExpanded === undefined) {
                window.venomzoneExpanded = false;
            }
            window.venomzoneExpanded = !window.venomzoneExpanded;
            extraRows.forEach(function(row) {
                row.style.display = window.venomzoneExpanded ? "table-row" : "none";
            });
            toggleRow.innerHTML = `<td colspan="2">${window.venomzoneExpanded ? "Show less..." : "Show more..."}</td>`;
        }

        // Cluster-specific classification
        function renderClusterClassification(nodeName) {
            // Get the cluster ID for the current node
            var clusterId = window.node_clusters ? window.node_clusters[nodeName] : null;
            console.log("Selected node:", nodeName, "has cluster id:", clusterId);
            if (clusterId === null) return;
            
            // Gather all nodes that belong to the same cluster
            var allNodes = Object.keys(window.node_clusters || {});
            var clusterNodes = allNodes.filter(function(n) {
                return window.node_clusters[n] === clusterId;
            });

            var perNodeAggregation = {};  // key: node, value: object { classification: meanE }
            
            clusterNodes.forEach(function(n) {
                var hits = window.vzone_hits && window.vzone_hits[n] ? window.vzone_hits[n] : [];
                var nodeMap = {};  // key: classification, value: array of evalues
                hits.forEach(function(item) {
                var cl = item.classification || "N/A";
                if (!nodeMap[cl]) {
                    nodeMap[cl] = [];
                }
                nodeMap[cl].push(item.evalue);
                });
                perNodeAggregation[n] = {};
                Object.keys(nodeMap).forEach(function(cl) {
                var values = nodeMap[cl];
                var sum = values.reduce(function(acc, val) { return acc + val; }, 0);
                var meanE = sum / values.length;
                perNodeAggregation[n][cl] = meanE;
                });
            });
            
            var aggregate = {};  // key: classification, value: { count: number, totalE: sum of per-node means }
            
            Object.keys(perNodeAggregation).forEach(function(n) {
                var nodeData = perNodeAggregation[n];
                Object.keys(nodeData).forEach(function(cl) {
                if (!aggregate[cl]) {
                    aggregate[cl] = { count: 0, totalE: 0 };
                }
                aggregate[cl].count += 1; // Count how many nodes have this classification
                aggregate[cl].totalE += nodeData[cl]; // Sum the mean E-values for this classification
                });
            });
            
            // Create rows for each classification
            var totalNodesInCluster = clusterNodes.length;
            var clusterRows = "";
            Object.keys(aggregate).forEach(function(key) {
                var count = aggregate[key].count;
                var avgE = (aggregate[key].totalE / count).toExponential(2);
                clusterRows += "<tr><td>" + key + "</td><td>" + count + " / " + totalNodesInCluster + "</td><td>" + avgE + "</td></tr>";
            });
            document.querySelector("#cluster-table tbody").innerHTML = clusterRows;
        }

        // Somewhat similar to the above, but for conotoxins
        function renderConotoxinNodeClassification(nodeName) {
            const matches = window.cono_hits[nodeName] || [];
            matches.sort((a,b)=>a.evalue - b.evalue);
            const tbody = document.querySelector("#conotoxin-node-table tbody");
            tbody.innerHTML = "";
            matches.slice(0,5).forEach(m=>{
                const tr = `<tr>
                            <td>${m.classification}</td>
                            <td>${m.evalue}</td>
                            </tr>`;
                tbody.innerHTML += tr;
            });
        }

        function renderConotoxinClusterClassification(nodeName) {
        const clusterId     = window.node_clusters[nodeName];
        const nodesInCluster = Object.keys(window.node_clusters)
                                    .filter(n => window.node_clusters[n] === clusterId);
        const perClassStats = {};  // { class: { count: x, totalE: y } }

        // For each node group hits by classification
        nodesInCluster.forEach(n => {
            const hits = window.cono_hits[n] || [];
            if (!hits.length) return;

            const byClass = hits.reduce((m, h) => {
            (m[h.classification] = m[h.classification] || []).push(h.evalue);
            return m;
            }, {});

            // For each class seen on this node, count votes
            Object.entries(byClass).forEach(([cls, eVals]) => {
            const meanE = eVals.reduce((s, v) => s + v, 0) / eVals.length;
            if (!perClassStats[cls]) perClassStats[cls] = { count: 0, totalE: 0 };
            perClassStats[cls].count  += 1;     // only one vote per node
            perClassStats[cls].totalE += meanE;
            });
        });

        const tbody       = document.querySelector("#conotoxin-cluster-table tbody");
        const totalNodes  = nodesInCluster.length;
        tbody.innerHTML   = "";

        Object.entries(perClassStats)
                .sort((a, b) => b[1].count - a[1].count)  // optional: sort by popularity
                .forEach(([cls, stats]) => {
            const avgE = (stats.totalE / stats.count).toExponential(2);
            tbody.innerHTML += `
            <tr>
                <td>${cls}</td>
                <td>${stats.count} / ${totalNodes}</td>
                <td>${avgE}</td>
            </tr>`;
        });
        }

        // Grab elements
        const dropZone    = document.getElementById("drop-zone");
        const fileInput   = document.getElementById("file-input");
        const infoDisplay = document.getElementById("dropped-file-info");
        const classifyBtn = document.getElementById("classify-btn");

        // When the drop zone is clicked, open the file dialog
        dropZone.addEventListener("click", () => fileInput.click());

        // Prevent default behavior for drag/drop
        ["dragenter","dragover","dragleave","drop"].forEach(evt => {
        dropZone.addEventListener(evt, e => {
            e.preventDefault();
            e.stopPropagation();
        });
        });

        // Highlight on dragover
        dropZone.addEventListener("dragover", () => {
        dropZone.classList.add("dragover");
        });
        dropZone.addEventListener("dragleave", () => {
        dropZone.classList.remove("dragover");
        });

        // Handle dropped files
        dropZone.addEventListener("drop", e => {
        dropZone.classList.remove("dragover");
        if (!e.dataTransfer.files.length) return;
        fileInput.files = e.dataTransfer.files;
        const f = e.dataTransfer.files[0];
        infoDisplay.textContent = `Selected file: ${f.name}`;
        classifyBtn.disabled = false;
        });

        // Also handle manual selection via dialog
        fileInput.addEventListener("change", () => {
        if (!fileInput.files.length) {
            infoDisplay.textContent = "No file selected.";
            classifyBtn.disabled = true;
            return;
        }
        const f = fileInput.files[0];
        infoDisplay.textContent = `Selected file: ${f.name}`;
        classifyBtn.disabled = false;
        });

        // Export information to CSV-file
        document.getElementById("export-btn").addEventListener("click", () => {
        const inverseSpeciesMapping = Object.entries(speciesMapping).reduce((inv, [latin, testId]) => {
            inv[testId] = latin;
            return inv;
        }, {});

        const nodeClusters = window.node_clusters;  // { nodeID: clusterID, ... }
        const clusters     = {};
        Object.entries(nodeClusters).forEach(([node, cid]) => {
            (clusters[cid] = clusters[cid] || []).push(node);
        });

        // Prepare CSV header and rows
        const header = [
            "cluster_id",
            "num_nodes",
            "species_representatives",
            "node_identifiers",
            "venomzone_cluster_hits",
            "conotoxin_cluster_hits",
            "comments"
        ].join(",");

        const lines = [ header ];

        Object.entries(clusters).forEach(([cid, nodes]) => {
            const numNodes = nodes.length;
            const latinSet = new Set(
            nodes.map(n => {
                const rec    = window.network.body.data.nodes.get(n);
                const testId = rec && rec.species;
                return inverseSpeciesMapping[testId] || testId;
            }).filter(Boolean)
            );
            let speciesReps = Array.from(latinSet).join(", ");
            if (latinSet.size > 1) {
            speciesReps = `(${speciesReps})`;
            }

            const nodeList = nodes.join(",");

            const vzAgg = {};
            nodes.forEach(n => {
            const hits = window.vzone_hits[n] || [];
            hits.forEach(h => {
                const cl = h.classification || "N/A";
                (vzAgg[cl] = vzAgg[cl] || { total: 0, count: 0 });
                vzAgg[cl].total += h.evalue;
                vzAgg[cl].count += 1;
            });
            });
            const vzStr = (() => {
                if (Object.keys(vzAgg).length === 0) return "N/A";
                return Object.entries(vzAgg)
                .map(([cl, { count, total }]) => {
                    const avg = (total / count).toExponential(2);
                    return `${cl}:${count}/${numNodes}:${avg}`;
                })
                .join("|");
            })();

            const ctAgg = {};
            nodes.forEach(n => {
            const hits = window.cono_hits[n] || [];
            const byClass = hits.reduce((m, h) => {
                (m[h.classification] = m[h.classification] || []).push(h.evalue);
                return m;
            }, {});
            Object.entries(byClass).forEach(([cl, vals]) => {
                const meanE = vals.reduce((s, v) => s + v, 0) / vals.length;
                (ctAgg[cl] = ctAgg[cl] || { total: 0, count: 0 });
                ctAgg[cl].total += meanE;
                ctAgg[cl].count += 1;
            });
            });
            const ctStr = (() => {
                if (Object.keys(ctAgg).length === 0) return "N/A";
                return Object.entries(ctAgg)
                .map(([cl, { count, total }]) => {
                    const avg = (total / count).toExponential(2);
                    return `${cl}:${count}/${numNodes}:${avg}`;
                })
                .join("|");
            })();

            // Comments always empty
            const comments = "";

            lines.push([
            cid,
            numNodes,
            `"${speciesReps}"`,
            `"${nodeList}"`,
            `"${vzStr}"`,
            `"${ctStr}"`,
            `"${comments}"`
            ].join(", "));
        });

        // Invoke download
        const blob = new Blob([ lines.join("\\n") ], { type: "text/csv" });
        const url  = URL.createObjectURL(blob);
        const a    = document.createElement("a");
        a.href     = url;
        a.download = "clusters_export.csv"; // Name of the downloaded file
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
        });

        // Zoom to a specific cluster upon clicking a cluster link or retrieving match from search
        function zoomToCluster(clusterId) {
        const rawNodeIds = Object.keys(window.node_clusters)
            .filter(nodeId => window.node_clusters[nodeId] == clusterId); // == to allow string/number
        const nodeIds = rawNodeIds.map(id => id);
        const validNodeIds = nodeIds.filter(id =>
            window.network.body.data.nodes.get(id) !== undefined
        );

        if (!validNodeIds.length) {
            console.warn(`No nodes found for cluster ${clusterId}`);
            return;
        }

        window.network.fit({
            nodes: validNodeIds,
            animation: {
            duration: 500,
            easingFunction: 'easeInOutQuad'
            }
        });
        }
    
        // On click Classify
        classifyBtn.addEventListener("click", () => {
            const f = fileInput.files[0];
            if (!f) {
                return alert("Please drop or select a file first.");
            }

            const fd = new FormData();
            fd.append("file", f);

            fetch("/classify", { method: "POST", body: fd })
            .then(r => r.json())
            .then(info => {
                const resultDiv = document.getElementById("classify-result");
                const hitsDiv   = document.getElementById("classify-summary");

                if (info.target) {
                resultDiv.innerHTML = `
                    <strong>Best match:</strong>
                    Cluster ${info.cluster}
                    | ${info.target}
                    (E-value = ${info.evalue}, bits = ${info.bitscore})
                `;
                } else {
                resultDiv.innerText = "No best match could be identified.";
                }

                if (info.hits && info.hits.length) {
                let html = `<h5>Top ${info.hits.length} hits:</h5>`;
                html += `<ol style="margin-left:1em">`;
                info.hits.forEach(h => {
                    html += `
                    <li>
                        <i>${h.target}</i>
                        | E-value = ${h.evalue},
                        bits = ${h.bitscore},
                        cluster = ${h.cluster ?? 'N/A'}
                    </li>`;
                });
                html += `</ol>`;
                hitsDiv.innerHTML = html;
                } else {
                hitsDiv.innerText = "No additional hits.";
                }

                // Reset all node colors
                const allNodes = window.network.body.data.nodes.get();
                const updates = allNodes.map(n => ({
                    id: n.id,
                    color: speciesColors[n.species] || n.color
                }));
                window.network.body.data.nodes.update(updates);

                // Highlight the returned cluster (if available)
                if (info.cluster !== undefined) {
                    const ids = Object.keys(window.node_clusters)
                                    .filter(n => window.node_clusters[n] === info.cluster);
                    ids.forEach(id => {
                        window.network.body.data.nodes.update({ id, color: "#000000" });
                    });
                zoomToCluster(info.cluster);
                }
            })
            .catch(err => {
                alert("Error: " + err);
            });
        });

        window.afdbDetailsCache = {};
        function getAFDBDetails(acc) {
        if (window.afdbDetailsCache[acc]) {
            return Promise.resolve(window.afdbDetailsCache[acc]);
        }
        return fetch(`https://rest.uniprot.org/uniprotkb/${acc}.json`)
            .then(r => r.json())
            .then(data => {
                const name = data.proteinDescription?.recommendedName?.fullName?.value || acc;
                let mf = "N/A"
                if (Array.isArray(data.keywords)) {
                    const mfKeywords = data.keywords.filter(kw => kw.category === "Molecular function");
                    if (mfKeywords.length) {
                        mf = mfKeywords.map(kw => kw.name).join(", ");
                    }
                }
                const details = { name, mf };
                window.afdbDetailsCache[acc] = details;
                return details;
            })
            .catch(() => {
            const details = { name: acc, mf: "N/A" };
            window.afdbDetailsCache[acc] = details;
            return details;
            });
        }

        // Buttons and input box for search functionality
        const btn = document.getElementById("search-btn");
        const box = document.getElementById("search-box");
        const results = document.getElementById("search-results");
        btn.addEventListener("click", () => { // Event listener for button click
        results.style.display = results.style.display === "block" ? "none" : "block";
        });

        // Handle enter or button click
        box.addEventListener("keypress", e => {
        if (e.key === "Enter") doSearch();
        });
        btn.addEventListener("click", doSearch);

        function doSearch() {
            const q = box.value.trim();
            if (!q) return clearResults();

            const clusterMatch = q.match(/^cluster\s*:\s*(\d+)$/i);
            let matchedNodes = [];

            if (clusterMatch) {
                // Cluster specific search
                const cid = clusterMatch[1];
                matchedNodes = window.network.body.data.nodes.get().filter(n =>
                String(window.node_clusters[n.id]) === cid
                ).map(n => ({
                id: n.id,
                term: `Cluster ${cid}`,
                cluster: cid
            }));
            renderResults(matchedNodes);
            return;
            } 

            // Classification specific search
            const term     = q.toLowerCase();
            const allNodes = window.network.body.data.nodes.get();
            const matched  = [];

            allNodes.forEach(n => {
            (window.vzone_hits[n.id] || []).forEach(entry => {
                const cls = entry.classification;
                if (cls.toLowerCase().includes(term)) {
                    matched.push({
                    id:      n.id,
                    term:    cls,
                    source:  "Venomzone",
                    evalue:  entry.evalue,
                    cluster: window.node_clusters[n.id]
                    });
                }
            });
            });

            allNodes.forEach(n => {
            (window.cono_hits[n.id] || []).forEach(entry => {
                const cls = entry.classification;
                if (cls.toLowerCase().includes(term)) {
                matched.push({
                    id:      n.id,
                    term:    cls,
                    source:  "Conotoxin",
                    evalue:  entry.evalue,
                    cluster: window.node_clusters[n.id]
                });
                }
            });
            });

            const afPromises = [];
            Object.entries(window.afdb_hits || {}).forEach(([nodeId, entries]) => {
            entries.forEach(entry => {
                // extract accession from entry.target
                const acc = extractUniprotId(entry.target);
                afPromises.push(
                getAFDBDetails(acc).then(details => {
                    const text = (details.name + " " + details.mf).toLowerCase();
                    if (text.includes(term)) {
                    matched.push({
                        id:      nodeId,
                        term:    `${details.name} [GO: ${details.mf}]`,
                        source:  "AFDB Proteome",
                        evalue:  entry.evalue,
                        cluster: window.node_clusters[nodeId]
                    });
                    }
                })
                );
            });
            });

            Promise.all(afPromises).then(() => renderResults(matched));
            return;
        }

        // Clear the dropdown
        function clearResults() {
            results.innerHTML = "";
            results.style.display = "none";
        }

        // Render matched nodes list
        function renderResults(nodes) {
            results.innerHTML = ""; 
            if (nodes.length === 0) {
                results.innerHTML = `<div>No matches found</div>`;
                return;
        }
        nodes.forEach(n => {
            const row = document.createElement("div");
            row.className = "result-item";

            // Node ID
            const idSpan = document.createElement("span");
            idSpan.innerText = n.id;

            // Exact matched term
            const termSpan = document.createElement("span");
            termSpan.className = "matched-term";
            termSpan.innerText = n.term;

            // Go to
            const go = document.createElement("span");
            go.className = "go-to";
            go.innerText = "Go to";
            go.addEventListener("click", () => {
                zoomToCluster(n.cluster);
            });

            row.appendChild(idSpan);
            row.appendChild(termSpan);
            row.appendChild(go);
            results.appendChild(row);
            });
            results.style.display = "block";
        }

        // Extract cluster ID from a node ID and display it
        window.network.on("selectNode", params => {
        if (params.nodes.length > 0) {
            const nodeId = params.nodes[0];
            const clusterId = window.node_clusters[nodeId];
            document.getElementById("current-cluster").innerText = clusterId ?? "Unknown";
        }
        });

        // Clear when pressing on blank space
        window.network.on("deselectNode", () => {
        document.getElementById("current-cluster").innerText = "None";
        });