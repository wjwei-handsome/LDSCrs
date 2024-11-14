use anyhow::Result;
use std::{
    fs::File,
    io::{BufRead, BufReader, Read},
    path::Path,
};

const BUFFER_SIZE: usize = 32 * 1024;

const MAGIC_MAX_LEN: usize = 6;
// compressed file magic number, ref: https://docs.rs/infer/latest/infer/archive/index.html
const GZ_MAGIC: [u8; 3] = [0x1f, 0x8b, 0x08];
const BZ_MAGIC: [u8; 3] = [0x42, 0x5a, 0x68];
const XZ_MAGIC: [u8; 6] = [0xfd, 0x37, 0x7a, 0x58, 0x5A, 0x00];

fn get_magic_num(path: &str) -> Result<[u8; MAGIC_MAX_LEN]> {
    let mut buffer: [u8; MAGIC_MAX_LEN] = [0; MAGIC_MAX_LEN];
    let mut fp = File::open(path)?;
    let _ = fp.read(&mut buffer)?;
    Ok(buffer)
}

fn is_gzipped(path: &str) -> Result<bool> {
    let buffer = get_magic_num(path)?;
    let gz_or_not =
        buffer[0] == GZ_MAGIC[0] && buffer[1] == GZ_MAGIC[1] && buffer[2] == GZ_MAGIC[2];
    Ok(gz_or_not || Path::new(path).extension().is_some_and(|ext| ext == "gz"))
}

fn is_bzipped(path: &str) -> Result<bool> {
    let buffer = get_magic_num(path)?;
    let bz_or_not =
        buffer[0] == BZ_MAGIC[0] && buffer[1] == BZ_MAGIC[1] && buffer[2] == BZ_MAGIC[2];
    Ok(bz_or_not || Path::new(path).extension().is_some_and(|ext| ext == "bz2"))
}

fn is_xz(path: &str) -> Result<bool> {
    let buffer = get_magic_num(path)?;
    let xz_or_not = buffer[0] == XZ_MAGIC[0]
        && buffer[1] == XZ_MAGIC[1]
        && buffer[2] == XZ_MAGIC[2]
        && buffer[3] == XZ_MAGIC[3]
        && buffer[4] == XZ_MAGIC[4]
        && buffer[5] == XZ_MAGIC[5];
    Ok(xz_or_not || Path::new(path).extension().is_some_and(|ext| ext == "xz"))
}

pub fn get_input_reader(path: &str) -> Result<Box<dyn BufRead + Send>> {
    let reader: Box<dyn BufRead + Send> = match File::open(path) {
        Ok(file) => {
            if is_xz(path)? {
                // decode xz compressed file
                Box::new(BufReader::with_capacity(
                    BUFFER_SIZE,
                    xz2::read::XzDecoder::new_multi_decoder(file),
                ))
            } else if is_gzipped(path)? {
                // decode gzip compressed file
                Box::new(BufReader::with_capacity(
                    BUFFER_SIZE,
                    flate2::read::MultiGzDecoder::new(file),
                ))
            } else if is_bzipped(path)? {
                // decode bzip2 compressed file
                Box::new(BufReader::with_capacity(
                    BUFFER_SIZE,
                    bzip2::read::MultiBzDecoder::new(file),
                ))
            } else {
                // stdin flag "-" covered
                Box::new(BufReader::with_capacity(BUFFER_SIZE, file))
            }
        }
        Err(_) => {
            anyhow::bail!("File not found: {:?}", path);
        }
    };

    Ok(reader)
}
