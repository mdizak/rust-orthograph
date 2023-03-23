use std::path::Path;
use std::fs;
use uuid::Uuid;
use std::fs::File;
use std::io::Write;
use std::env::temp_dir;

pub fn create_dir(dirname: &String) {

    if !Path::new(&dirname).exists() {
        match fs::create_dir_all(&dirname) {
            Ok(_dir) => {},
            Err(_error) => panic!("Unable to create directory at {}.", dirname)
        };
    }

}

pub fn remove_dir(dirname: &String) {

    if Path::new(&dirname).exists() {
        match fs::remove_dir_all(&dirname) {
            Ok(_dir) => {},
            Err(_error) => panic!("Unable to remove directory at {}.", dirname)
        };
    }

}

pub fn recreate_dir(dirname: &String) {
    remove_dir(dirname);
    create_dir(dirname);
}

pub fn create_tmp_file(contents: &String) -> String {

    // Get filename
    let filename = gen_tmp_filename();
    let path = Path::new(&filename);

    // Open file
    let mut fh = match File::create(&path) {
        Ok(res) => res,
        Err(e) => panic!("Unable to open temporary file at {}, error: {}", filename, e)
    };

    // Write to file
    match fh.write_all(&contents.as_bytes()) {
        Ok(res) => res,
        Err(e) => panic!("Unable to write to temporary file {}, error: {}", filename, e)
    };

    filename
}

pub fn gen_tmp_filename() -> String {

    // Get filename
    let dir = temp_dir();
    let filename = format!("{}/{}", dir.to_str().unwrap(), Uuid::new_v4());
    //let filename = format!("/dev/shm/{}", Uuid::new_v4());

    filename
}

pub fn open_file(filename: String) -> File {

    let path = Path::new(&filename);
    let fh = match File::create(&path) {
        Ok(res) => res,
        Err(e) => panic!("Unable to open file for writing, {}, error: {}", filename, e)
    };

    fh
}


