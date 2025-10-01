# GitHub Actions Setup

This repository uses GitHub Actions to automatically build and push Docker images to DockerHub when certain files change.

## Workflows

### 1. Docker Test (`docker-test.yml`)
- Triggers on pull requests to `develop` and `main` branches when any of these files/directories change:
  - `Dockerfile`
  - `bin/` directory
  - `panel_information/` directory
- Builds the Docker image locally to test that it builds successfully
- Does NOT push to DockerHub (no credentials needed)

### 2. Docker Deploy (`docker-deploy.yml`)
- Triggers on pushes to `develop` and `main` branches when any of these files/directories change:
  - `Dockerfile`
  - `bin/` directory
  - `panel_information/` directory
- Builds and pushes the Docker image to DockerHub with appropriate tags:
  - `dev` tag for the develop branch
  - `latest` tag for the main branch
  - Branch-specific tags with commit SHA

## Image Tags

The workflow will create images with the following naming convention:
- `yourusername/mad4hatter:dev` - for develop branch
- `yourusername/mad4hatter:latest` - for main branch
- `yourusername/mad4hatter:develop-<commit-sha>` - for specific commits on develop
- `yourusername/mad4hatter:main-<commit-sha>` - for specific commits on main

## Usage

Once the secrets are configured, the workflow will automatically run whenever:
- You push changes to the develop or main branch that affect the Dockerfile, bin/, or panel_information/ directories (builds and pushes to DockerHub)
- A pull request is opened or updated that affects these files (builds locally to test, no push)

The workflow includes caching to speed up subsequent builds.
