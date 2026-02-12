<div class="title-logo-container">
  <div class="title-section">
    <h1>The MAD4HATTER Amplicon Sequencing Pipeline</h1>
    <p><strong>M</strong>ultiplex <strong>A</strong>mplicons for <strong>D</strong>rugs, <strong>D</strong>iagnostics, <strong>D</strong>iversity, and <strong>D</strong>ifferentiation using <strong>H</strong>igh <strong>T</strong>hroughput <strong>T</strong>argeted <strong>R</strong>esequencing</p>
  </div>
  <div class="logo-section">
    <img src="assets/logo.svg" alt="MAD4HATTER Logo" class="home-logo">
  </div>
</div>

Welcome to the MAD4HATTER pipeline!

Mad4hatter is a bioinformatics pipeline designed to analyse Plasmodium Illumina Amplicon sequencing data. It processes raw FASTQ files and produces an allele table, core QC metrics, and drug-resistance information. It was originally developed for the MAD4HatTeR panel but has since been adapted to support additional panels. While the pipeline can be run on any panel, it was optimised using MAD4HatTeR data; panels with substantially different properties, such as very short amplicon targets, may require additional tuning to achieve optimal performance. Several commonly used panels are preconfigured for convenience, and new panels can be easily added through simple configuration.

The pipeline can be run using **Nextflow** (command-line interface) or **Terra** (cloud-based platform with a graphical interface).

**Choose Your Platform**

<div class="platform-buttons">
  <a href="getting-started/" class="platform-button">
    <img src="assets/images/nextflow_logo.png" alt="Nextflow Logo" class="platform-logo">
    <span class="platform-text">Run with Nextflow</span>
    <span class="platform-description">Command-line interface for local computers and HPC clusters</span>
  </a>
  
  <a href="https://app.terra.bio/#workspaces/gates-malaria/Mad4Hatter" target="_blank" rel="noopener" class="platform-button">
    <img src="assets/images/terra-logo.png" alt="Terra Logo" class="platform-logo">
    <span class="platform-text">Run on Terra</span>
    <span class="platform-description">Cloud-based platform with graphical interface - no installation required</span>
  </a>
</div> 

**Thank you to the contributors of the MAD4HATTER pipeline!**

<ul class="contributors-list" id="contributors-list">
  <li>Loading contributors...</li>
</ul>

<script>
// Fetch GitHub contributors
fetch('https://api.github.com/repos/eppicenter/mad4hatter/contributors')
  .then(response => response.json())
  .then(contributors => {
    const list = document.getElementById('contributors-list');
    if (contributors && contributors.length > 0) {
      list.innerHTML = contributors.map(contributor => `
        <li>
          <a href="${contributor.html_url}" title="${contributor.login}" target="_blank" rel="noopener">
            <img src="${contributor.avatar_url}" alt="${contributor.login}" loading="lazy">
          </a>
        </li>
      `).join('');
    } else {
      list.innerHTML = '<li>No contributors found.</li>';
    }
  })
  .catch(error => {
    console.error('Error fetching contributors:', error);
    document.getElementById('contributors-list').innerHTML = '<li>Unable to load contributors.</li>';
  });
</script>
