use std::env;
use std::fs;
use std::process::Command;

fn set_git_version(){
    let version = env::var("CARGO_PKG_VERSION").unwrap();

    let child = Command::new("git").args(["describe", "--always"]).output();
    match child {
        Ok(child) => {
            let buf = String::from_utf8(child.stdout).expect("failed to read stdout");
            println!("cargo:rustc-env=VERSION={version}-{buf}");
        }
        Err(err) => {
            eprintln!("`git describe` err: {err}");
            println!("cargo:rustc-env=VERSION={version}");
        }
    }
}

fn set_update_checker() {
    let manifest_dir = env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    let git_config = format!("{}/.git/config", manifest_dir);

    let (mut owner, mut repo) = ("unknown".to_string(), "unknown".to_string());

    if let Ok(content) = fs::read_to_string(&git_config) {
        if let Some(line) = content.lines().find(|l| l.trim().starts_with("url =")) {
            let url = line.trim()["url =".len()..].trim();

            // 1. git@github.com:owner/repo.git
            // 2. https://github.com/owner/repo.git
            // 3. git://github.com/owner/repo.git
            let repo_part = if let Some(idx) = url.find("github.com") {
                &url[idx + "github.com/".len()..]
            }else {
                ""
            };

            let repo_part = repo_part.trim_end_matches(".git");

            let mut parts = repo_part.splitn(2, '/');
            if let (Some(o), Some(r)) = (parts.next(), parts.next()) {
                if !o.is_empty() && !r.is_empty() {
                    owner = o.to_string();
                    repo = r.to_string();
                }
            }
        }
    }
    println!("cargo:rustc-env=GIT_OWNER={owner}");
    println!("cargo:rustc-env=GIT_REPO={repo}");
}

fn main() {
    set_git_version();
    set_update_checker();
}
