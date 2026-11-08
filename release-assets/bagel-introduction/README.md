# Upload the Bagel introduction video

Status: prepared locally; the video has not been uploaded. The vignette's video link will work after the release below is published.

This folder is outside `bagelR/`, so it is not included in the CRAN package.

## Files

* `bagel-introduction.mp4`: the video to attach to the release.
* `release-notes.md`: text to paste into the release description, including credits.

## Publish on the existing GitHub repository

You need write access to `grosed/bagel`.

1. Open [Draft a new release](https://github.com/grosed/bagel/releases/new).
2. Choose or create the tag **`bagel-introduction`**, targeting the repository's default branch. Use this exact tag: the vignette link includes it.
3. Set the title to **Bagel introduction video** and paste the contents of `release-notes.md` into the description.
4. Attach **`bagel-introduction.mp4`** in the release assets area, keeping the filename unchanged. Wait for the upload to finish.
5. Leave **Set as latest release** unselected, since this release distributes a video rather than a new package version. Publish the release after the asset is attached. A draft is not publicly accessible.
6. Open the video link below while signed out and confirm it plays or downloads. Do this before submitting the R package to CRAN.

## Video address used by the vignette

https://github.com/grosed/bagel/releases/download/bagel-introduction/bagel-introduction.mp4

Keep the tag and filename above to preserve this address. If you use a different repository, tag or filename, update both links in `bagelR/vignettes/bagelR.Rmd` and rebuild the package.

GitHub's [release instructions](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository) describe creating a release and attaching files.
