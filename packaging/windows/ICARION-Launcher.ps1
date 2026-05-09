Add-Type -AssemblyName System.Windows.Forms
Add-Type -AssemblyName System.Drawing

$ErrorActionPreference = "Stop"
[System.Windows.Forms.Application]::EnableVisualStyles()
[System.Windows.Forms.Application]::SetCompatibleTextRenderingDefault($false)

$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$IcarionExe = Join-Path $Root "bin\icarion.exe"
$DefaultConfigDir = Join-Path $Root "examples"
$RunLogDir = Join-Path $Root "launcher-logs"
if (Test-Path (Join-Path $Root "analysis")) {
    $AnalysisDir = Join-Path $Root "analysis"
} else {
    $AnalysisDir = Join-Path $Root "share\icarion\analysis"
}

$ColorBackground = [System.Drawing.Color]::FromArgb(246, 248, 250)
$ColorHeader = [System.Drawing.Color]::FromArgb(32, 43, 54)
$ColorText = [System.Drawing.Color]::FromArgb(25, 32, 39)
$ColorSubtle = [System.Drawing.Color]::FromArgb(96, 105, 114)
$ColorBorder = [System.Drawing.Color]::FromArgb(210, 216, 222)
$ColorPrimary = [System.Drawing.Color]::FromArgb(33, 116, 184)
$ColorStop = [System.Drawing.Color]::FromArgb(178, 62, 62)
$ColorConsole = [System.Drawing.Color]::FromArgb(18, 24, 31)
$ColorConsoleText = [System.Drawing.Color]::FromArgb(226, 232, 240)
$FontUi = New-Object System.Drawing.Font("Segoe UI", 9)
$FontTitle = New-Object System.Drawing.Font("Segoe UI Semibold", 16)
$FontSmall = New-Object System.Drawing.Font("Segoe UI", 8.5)
$FontMono = New-Object System.Drawing.Font("Consolas", 9)

function Show-Message {
    param(
        [string]$Message,
        [string]$Title = "ICARION Launcher",
        [System.Windows.Forms.MessageBoxIcon]$Icon = [System.Windows.Forms.MessageBoxIcon]::Information
    )
    [System.Windows.Forms.MessageBox]::Show(
        $form,
        $Message,
        $Title,
        [System.Windows.Forms.MessageBoxButtons]::OK,
        $Icon
    ) | Out-Null
}

function Append-Log {
    param([string]$Text)
    if ([string]::IsNullOrWhiteSpace($Text)) {
        return
    }
    $logBox.AppendText($Text + [Environment]::NewLine)
    $logBox.SelectionStart = $logBox.Text.Length
    $logBox.ScrollToCaret()
}

function Read-LogText {
    param([string]$Path)

    if (-not (Test-Path -LiteralPath $Path)) {
        return ""
    }

    $stream = [System.IO.File]::Open($Path, [System.IO.FileMode]::Open, [System.IO.FileAccess]::Read, [System.IO.FileShare]::ReadWrite)
    try {
        $reader = New-Object System.IO.StreamReader($stream, [System.Text.Encoding]::UTF8, $true)
        try {
            $text = $reader.ReadToEnd()
        } finally {
            $reader.Dispose()
        }
    } finally {
        $stream.Dispose()
    }
    $escape = [regex]::Escape([string][char]27)
    return [regex]::Replace($text, "$escape\[[0-9;?]*[ -/]*[@-~]", "")
}

function Resolve-ConfigPath {
    param([string]$PathText)

    $candidate = $PathText.Trim().Trim('"')
    if ([string]::IsNullOrWhiteSpace($candidate)) {
        throw "Please select a JSON config file."
    }
    if (-not (Test-Path -LiteralPath $candidate)) {
        throw "This path does not exist:`n$candidate"
    }
    $item = Get-Item -LiteralPath $candidate
    if ($item.PSIsContainer) {
        throw "You selected a folder, not a config file.`n`nOpen the folder and select a .json config, for example examples\ims\ims_basic.json."
    }
    if ($item.Extension -ne ".json") {
        throw "ICARION configs must be .json files.`n`nSelected:`n$candidate"
    }
    return $item.FullName
}

function Get-PackageRelativePath {
    param([string]$PathText)

    $fullPath = [System.IO.Path]::GetFullPath($PathText)
    $rootPath = [System.IO.Path]::GetFullPath($Root).TrimEnd('\', '/')
    $prefix = $rootPath + [System.IO.Path]::DirectorySeparatorChar
    if ($fullPath.StartsWith($prefix, [System.StringComparison]::OrdinalIgnoreCase)) {
        return $fullPath.Substring($prefix.Length)
    }
    return $fullPath
}

function Select-TrajectoryFile {
    $dialog = New-Object System.Windows.Forms.OpenFileDialog
    $dialog.Title = "Select ICARION trajectory file"
    $dialog.Filter = "HDF5 trajectory (*.h5;*.hdf5)|*.h5;*.hdf5|All files (*.*)|*.*"
    $dialog.CheckFileExists = $true
    $dialog.CheckPathExists = $true
    $dialog.InitialDirectory = $Root
    if ($dialog.ShowDialog($form) -eq [System.Windows.Forms.DialogResult]::OK) {
        return $dialog.FileName
    }
    return $null
}

function Show-ImagePreview {
    param(
        [string]$ImagePath,
        [string]$Title = "ICARION Analysis"
    )
    if (-not (Test-Path -LiteralPath $ImagePath)) {
        return
    }

    $preview = New-Object System.Windows.Forms.Form
    $preview.Text = $Title
    $preview.StartPosition = "CenterParent"
    $preview.Size = New-Object System.Drawing.Size(1000, 720)
    $preview.MinimumSize = New-Object System.Drawing.Size(720, 480)
    $preview.BackColor = $ColorBackground
    $preview.KeyPreview = $true

    $buttonPanel = New-Object System.Windows.Forms.Panel
    $buttonPanel.Dock = [System.Windows.Forms.DockStyle]::Bottom
    $buttonPanel.Height = 44
    $buttonPanel.Padding = New-Object System.Windows.Forms.Padding(8)
    $buttonPanel.BackColor = [System.Drawing.Color]::White

    $closePreviewButton = New-FlatButton "Close" $ColorPrimary
    $closePreviewButton.Dock = [System.Windows.Forms.DockStyle]::Right
    $closePreviewButton.Width = 86
    $closePreviewButton.Add_Click({
        param($sender, $eventArgs)
        $sender.FindForm().Close()
    })

    $openImageButton = New-FlatButton "Open" ([System.Drawing.Color]::FromArgb(229, 234, 240)) $ColorText
    $openImageButton.Dock = [System.Windows.Forms.DockStyle]::Right
    $openImageButton.Width = 86
    $openImageButton.Add_Click({ Invoke-Item -LiteralPath $ImagePath })

    $openFolderButton = New-FlatButton "Folder" ([System.Drawing.Color]::FromArgb(229, 234, 240)) $ColorText
    $openFolderButton.Dock = [System.Windows.Forms.DockStyle]::Right
    $openFolderButton.Width = 86
    $openFolderButton.Add_Click({ Invoke-Item -LiteralPath (Split-Path -Parent $ImagePath) })

    $buttonPanel.Controls.Add($closePreviewButton)
    $buttonPanel.Controls.Add($openImageButton)
    $buttonPanel.Controls.Add($openFolderButton)

    $picture = New-Object System.Windows.Forms.PictureBox
    $picture.Dock = [System.Windows.Forms.DockStyle]::Fill
    $picture.SizeMode = [System.Windows.Forms.PictureBoxSizeMode]::Zoom
    $picture.BackColor = [System.Drawing.Color]::White
    $picture.Load($ImagePath)

    $pathLabel = New-Object System.Windows.Forms.TextBox
    $pathLabel.Dock = [System.Windows.Forms.DockStyle]::Bottom
    $pathLabel.Height = 28
    $pathLabel.ReadOnly = $true
    $pathLabel.BorderStyle = [System.Windows.Forms.BorderStyle]::FixedSingle
    $pathLabel.Text = $ImagePath

    $preview.Controls.Add($picture)
    $preview.Controls.Add($pathLabel)
    $preview.Controls.Add($buttonPanel)
    $preview.AcceptButton = $closePreviewButton
    $preview.Add_KeyDown({
        param($sender, $eventArgs)
        if ($eventArgs.KeyCode -eq [System.Windows.Forms.Keys]::Escape) {
            $sender.Close()
        }
    })
    [void]$preview.Show($form)
}

function Get-PythonCommand {
    if (Get-Command py.exe -ErrorAction SilentlyContinue) {
        return "py -3"
    }
    if (Get-Command python.exe -ErrorAction SilentlyContinue) {
        return "python"
    }
    if (Get-Command python3.exe -ErrorAction SilentlyContinue) {
        return "python3"
    }
    return $null
}

function Set-RunControlsEnabled {
    param([bool]$Enabled)
    $runButton.Enabled = $Enabled
    $browseButton.Enabled = $Enabled
    $examplesButton.Enabled = $Enabled
    $analysisMobilityButton.Enabled = $Enabled
    $analysisArrivalButton.Enabled = $Enabled
    $analysisTrajectoryButton.Enabled = $Enabled
    $analysisMeanButton.Enabled = $Enabled
    $analysisEliminationButton.Enabled = $Enabled
    $analysisAnimateButton.Enabled = $Enabled
    $stopButton.Enabled = -not $Enabled
}

function New-FlatButton {
    param(
        [string]$Text,
        [System.Drawing.Color]$BackColor,
        [System.Drawing.Color]$ForeColor = [System.Drawing.Color]::White
    )
    $button = New-Object System.Windows.Forms.Button
    $button.Text = $Text
    $button.FlatStyle = [System.Windows.Forms.FlatStyle]::Flat
    $button.FlatAppearance.BorderSize = 0
    $button.BackColor = $BackColor
    $button.ForeColor = $ForeColor
    $button.Font = $FontUi
    $button.Cursor = [System.Windows.Forms.Cursors]::Hand
    return $button
}

$form = New-Object System.Windows.Forms.Form
$form.Text = "ICARION Launcher"
$form.StartPosition = "CenterScreen"
$form.Size = New-Object System.Drawing.Size(900, 620)
$form.MinimumSize = New-Object System.Drawing.Size(900, 500)
$form.BackColor = $ColorBackground
$form.Font = $FontUi

$headerPanel = New-Object System.Windows.Forms.Panel
$headerPanel.Dock = [System.Windows.Forms.DockStyle]::Top
$headerPanel.Height = 86
$headerPanel.BackColor = $ColorHeader

$titleLabel = New-Object System.Windows.Forms.Label
$titleLabel.Text = "ICARION"
$titleLabel.ForeColor = [System.Drawing.Color]::White
$titleLabel.Font = $FontTitle
$titleLabel.Location = New-Object System.Drawing.Point(20, 16)
$titleLabel.AutoSize = $true

$subtitleLabel = New-Object System.Windows.Forms.Label
$subtitleLabel.Text = "Load a JSON config and run an ion simulation"
$subtitleLabel.ForeColor = [System.Drawing.Color]::FromArgb(207, 216, 225)
$subtitleLabel.Font = $FontUi
$subtitleLabel.Location = New-Object System.Drawing.Point(22, 50)
$subtitleLabel.AutoSize = $true

$headerPanel.Controls.AddRange(@($titleLabel, $subtitleLabel))

$contentPanel = New-Object System.Windows.Forms.Panel
$contentPanel.Dock = [System.Windows.Forms.DockStyle]::Fill
$contentPanel.Padding = New-Object System.Windows.Forms.Padding(18, 18, 18, 14)
$contentPanel.BackColor = $ColorBackground

$configPanel = New-Object System.Windows.Forms.Panel
$configPanel.Dock = [System.Windows.Forms.DockStyle]::Top
$configPanel.Height = 202
$configPanel.BackColor = [System.Drawing.Color]::White
$configPanel.Padding = New-Object System.Windows.Forms.Padding(16, 12, 16, 12)

$configPanel.add_Paint({
    param($sender, $eventArgs)
    $rect = New-Object System.Drawing.Rectangle(0, 0, ($sender.Width - 1), ($sender.Height - 1))
    $pen = New-Object System.Drawing.Pen($ColorBorder)
    $eventArgs.Graphics.DrawRectangle($pen, $rect)
    $pen.Dispose()
})

$configLabel = New-Object System.Windows.Forms.Label
$configLabel.Text = "Config file"
$configLabel.ForeColor = $ColorText
$configLabel.Font = New-Object System.Drawing.Font("Segoe UI Semibold", 9)
$configLabel.Location = New-Object System.Drawing.Point(16, 12)
$configLabel.AutoSize = $true

$configBox = New-Object System.Windows.Forms.TextBox
$configBox.Location = New-Object System.Drawing.Point(16, 40)
$configBox.Size = New-Object System.Drawing.Size(820, 24)
$configBox.Anchor = [System.Windows.Forms.AnchorStyles]::Top -bor [System.Windows.Forms.AnchorStyles]::Left
$configBox.Font = $FontUi
$configBox.BorderStyle = [System.Windows.Forms.BorderStyle]::FixedSingle

$browseButton = New-FlatButton "Browse..." ([System.Drawing.Color]::FromArgb(229, 234, 240)) $ColorText
$browseButton.Location = New-Object System.Drawing.Point(16, 78)
$browseButton.Size = New-Object System.Drawing.Size(86, 30)
$browseButton.Anchor = [System.Windows.Forms.AnchorStyles]::Top -bor [System.Windows.Forms.AnchorStyles]::Left

$examplesButton = New-FlatButton "Examples..." ([System.Drawing.Color]::FromArgb(229, 234, 240)) $ColorText
$examplesButton.Location = New-Object System.Drawing.Point(110, 78)
$examplesButton.Size = New-Object System.Drawing.Size(96, 30)
$examplesButton.Anchor = [System.Windows.Forms.AnchorStyles]::Top -bor [System.Windows.Forms.AnchorStyles]::Left

$runButton = New-FlatButton "Run" $ColorPrimary
$runButton.Location = New-Object System.Drawing.Point(214, 78)
$runButton.Size = New-Object System.Drawing.Size(68, 30)
$runButton.Anchor = [System.Windows.Forms.AnchorStyles]::Top -bor [System.Windows.Forms.AnchorStyles]::Left

$stopButton = New-FlatButton "Stop" $ColorStop
$stopButton.Text = "Stop"
$stopButton.Enabled = $false
$stopButton.Location = New-Object System.Drawing.Point(290, 78)
$stopButton.Size = New-Object System.Drawing.Size(60, 30)
$stopButton.Anchor = [System.Windows.Forms.AnchorStyles]::Top -bor [System.Windows.Forms.AnchorStyles]::Left

$analysisMobilityButton = New-FlatButton "IMS mobility..." ([System.Drawing.Color]::FromArgb(61, 133, 108))
$analysisMobilityButton.Location = New-Object System.Drawing.Point(16, 116)
$analysisMobilityButton.Size = New-Object System.Drawing.Size(116, 30)
$analysisMobilityButton.Anchor = [System.Windows.Forms.AnchorStyles]::Top -bor [System.Windows.Forms.AnchorStyles]::Left

$analysisArrivalButton = New-FlatButton "Arrival times..." ([System.Drawing.Color]::FromArgb(92, 105, 130))
$analysisArrivalButton.Location = New-Object System.Drawing.Point(140, 116)
$analysisArrivalButton.Size = New-Object System.Drawing.Size(118, 30)
$analysisArrivalButton.Anchor = [System.Windows.Forms.AnchorStyles]::Top -bor [System.Windows.Forms.AnchorStyles]::Left

$analysisTrajectoryButton = New-FlatButton "Trajectories..." ([System.Drawing.Color]::FromArgb(92, 105, 130))
$analysisTrajectoryButton.Location = New-Object System.Drawing.Point(266, 116)
$analysisTrajectoryButton.Size = New-Object System.Drawing.Size(112, 30)
$analysisTrajectoryButton.Anchor = [System.Windows.Forms.AnchorStyles]::Top -bor [System.Windows.Forms.AnchorStyles]::Left

$analysisMeanButton = New-FlatButton "Mean positions..." ([System.Drawing.Color]::FromArgb(92, 105, 130))
$analysisMeanButton.Location = New-Object System.Drawing.Point(386, 116)
$analysisMeanButton.Size = New-Object System.Drawing.Size(124, 30)
$analysisMeanButton.Anchor = [System.Windows.Forms.AnchorStyles]::Top -bor [System.Windows.Forms.AnchorStyles]::Left

$analysisEliminationButton = New-FlatButton "Eliminations..." ([System.Drawing.Color]::FromArgb(92, 105, 130))
$analysisEliminationButton.Location = New-Object System.Drawing.Point(518, 116)
$analysisEliminationButton.Size = New-Object System.Drawing.Size(116, 30)
$analysisEliminationButton.Anchor = [System.Windows.Forms.AnchorStyles]::Top -bor [System.Windows.Forms.AnchorStyles]::Left

$analysisAnimateButton = New-FlatButton "Animate..." ([System.Drawing.Color]::FromArgb(92, 105, 130))
$analysisAnimateButton.Location = New-Object System.Drawing.Point(642, 116)
$analysisAnimateButton.Size = New-Object System.Drawing.Size(96, 30)
$analysisAnimateButton.Anchor = [System.Windows.Forms.AnchorStyles]::Top -bor [System.Windows.Forms.AnchorStyles]::Left

$analysisHintLabel = New-Object System.Windows.Forms.Label
$analysisHintLabel.Text = "Analysis uses Python packages from requirements-analysis.txt"
$analysisHintLabel.ForeColor = $ColorSubtle
$analysisHintLabel.Font = $FontSmall
$analysisHintLabel.Location = New-Object System.Drawing.Point(16, 158)
$analysisHintLabel.AutoSize = $true

$configPanel.Controls.AddRange(@(
    $configLabel,
    $configBox,
    $browseButton,
    $examplesButton,
    $runButton,
    $stopButton,
    $analysisMobilityButton,
    $analysisArrivalButton,
    $analysisTrajectoryButton,
    $analysisMeanButton,
    $analysisEliminationButton,
    $analysisAnimateButton,
    $analysisHintLabel
))

$statusPanel = New-Object System.Windows.Forms.Panel
$statusPanel.Dock = [System.Windows.Forms.DockStyle]::Bottom
$statusPanel.Height = 32
$statusPanel.BackColor = [System.Drawing.Color]::White

$statusPanel.add_Paint({
    param($sender, $eventArgs)
    $pen = New-Object System.Drawing.Pen($ColorBorder)
    $eventArgs.Graphics.DrawLine($pen, 0, 0, $sender.Width, 0)
    $pen.Dispose()
})

$statusLabel = New-Object System.Windows.Forms.Label
$statusLabel.Text = "Ready"
$statusLabel.ForeColor = $ColorSubtle
$statusLabel.Font = $FontSmall
$statusLabel.Dock = [System.Windows.Forms.DockStyle]::Fill
$statusLabel.Padding = New-Object System.Windows.Forms.Padding(12, 8, 12, 0)
$statusPanel.Controls.Add($statusLabel)

$logBox = New-Object System.Windows.Forms.TextBox
$logBox.Dock = [System.Windows.Forms.DockStyle]::Fill
$logBox.Margin = New-Object System.Windows.Forms.Padding(0, 14, 0, 10)
$logBox.Multiline = $true
$logBox.ScrollBars = "Both"
$logBox.WordWrap = $false
$logBox.ReadOnly = $true
$logBox.Font = $FontMono
$logBox.BorderStyle = [System.Windows.Forms.BorderStyle]::None
$logBox.BackColor = $ColorConsole
$logBox.ForeColor = $ColorConsoleText

$logPanel = New-Object System.Windows.Forms.Panel
$logPanel.Dock = [System.Windows.Forms.DockStyle]::Fill
$logPanel.Padding = New-Object System.Windows.Forms.Padding(12)
$logPanel.BackColor = $ColorConsole
$logPanel.Controls.Add($logBox)

$contentPanel.Controls.Add($configPanel)
$contentPanel.Controls.Add($logPanel)
$configPanel.BringToFront()

$form.Controls.Add($contentPanel)
$form.Controls.Add($statusPanel)
$form.Controls.Add($headerPanel)

$script:process = $null
$script:runLogPath = $null
$script:currentPlotPath = $null
$script:currentCsvPath = $null
$script:currentRunName = $null

$logTimer = New-Object System.Windows.Forms.Timer
$logTimer.Interval = 500
$logTimer.Add_Tick({
    try {
        if ($script:runLogPath -and (Test-Path -LiteralPath $script:runLogPath)) {
            $content = Read-LogText $script:runLogPath
            if ($null -ne $content -and $content -ne $logBox.Text) {
                $logBox.Text = $content
                $logBox.SelectionStart = $logBox.Text.Length
                $logBox.ScrollToCaret()
            }
        }

        if ($script:process -and $script:process.HasExited) {
            $code = $script:process.ExitCode
            $logTimer.Stop()
            Set-RunControlsEnabled $true
            $statusLabel.Text = "Finished with exit code $code"
            Append-Log ""
            Append-Log "Finished with exit code $code"
            if ($code -eq 0 -and $script:currentPlotPath -and (Test-Path -LiteralPath $script:currentPlotPath)) {
                Show-Message `
                    -Message "Analysis finished successfully.`n`nPlot:`n$script:currentPlotPath`n`nCSV:`n$script:currentCsvPath" `
                    -Title "Analysis finished" `
                    -Icon ([System.Windows.Forms.MessageBoxIcon]::Information)
                Show-ImagePreview $script:currentPlotPath $script:currentRunName
            }
            if ($code -ne 0) {
                Show-Message `
                    -Message "ICARION stopped with exit code $code.`n`nThe log was saved here:`n$script:runLogPath" `
                    -Title "Run failed" `
                    -Icon ([System.Windows.Forms.MessageBoxIcon]::Error)
            }
            $script:process.Dispose()
            $script:process = $null
        }
    } catch {
        $logTimer.Stop()
        Set-RunControlsEnabled $true
        $statusLabel.Text = "Launcher error"
        Append-Log $_.Exception.Message
        Show-Message `
            -Message $_.Exception.Message `
            -Title "Launcher error" `
            -Icon ([System.Windows.Forms.MessageBoxIcon]::Error)
    }
})

$browseButton.Add_Click({
    $dialog = New-Object System.Windows.Forms.OpenFileDialog
    $dialog.Title = "Select ICARION config"
    $dialog.Filter = "JSON config (*.json)|*.json"
    $dialog.CheckFileExists = $true
    $dialog.CheckPathExists = $true
    if (Test-Path $DefaultConfigDir) {
        $dialog.InitialDirectory = $DefaultConfigDir
    }
    if ($dialog.ShowDialog($form) -eq [System.Windows.Forms.DialogResult]::OK) {
        $configBox.Text = $dialog.FileName
    }
})

$examplesButton.Add_Click({
    $dialog = New-Object System.Windows.Forms.OpenFileDialog
    $dialog.Title = "Select an ICARION example config"
    $dialog.Filter = "JSON config (*.json)|*.json"
    $dialog.CheckFileExists = $true
    $dialog.CheckPathExists = $true
    if (Test-Path $DefaultConfigDir) {
        $dialog.InitialDirectory = $DefaultConfigDir
    }
    if ($dialog.ShowDialog($form) -eq [System.Windows.Forms.DialogResult]::OK) {
        $configBox.Text = $dialog.FileName
    }
})

$runButton.Add_Click({
    if (-not (Test-Path $IcarionExe)) {
        Show-Message `
            -Message "Cannot find bin\icarion.exe next to this launcher." `
            -Title "ICARION Launcher" `
            -Icon ([System.Windows.Forms.MessageBoxIcon]::Error)
        return
    }

    try {
        $configPath = Resolve-ConfigPath $configBox.Text
    } catch {
        Show-Message `
            -Message $_.Exception.Message `
            -Title "Invalid config selection" `
            -Icon ([System.Windows.Forms.MessageBoxIcon]::Warning)
        return
    }

    New-Item -ItemType Directory -Force $RunLogDir | Out-Null
    $runStamp = Get-Date -Format "yyyyMMdd-HHmmss"
    $runLogPath = Join-Path $RunLogDir "icarion-run-$runStamp.log"
    $runCmdPath = Join-Path $RunLogDir "icarion-run-$runStamp.cmd"
    $configArg = Get-PackageRelativePath $configPath
    $script:runLogPath = $runLogPath
    $script:currentPlotPath = $null
    $script:currentCsvPath = $null
    $script:currentRunName = $null

    $logBox.Clear()
    Append-Log "ICARION: $IcarionExe"
    Append-Log "Config:  $configPath"
    Append-Log "Argument: $configArg"
    Append-Log "Working directory: $Root"
    Append-Log "Wrapper: $runCmdPath"
    Append-Log "Log:     $runLogPath"
    Append-Log ""
    $utf8NoBom = New-Object System.Text.UTF8Encoding -ArgumentList $false
    [System.IO.File]::WriteAllText(
        $runLogPath,
        "ICARION: $IcarionExe`r`nConfig:  $configPath`r`nArgument: $configArg`r`nWorking directory: $Root`r`nWrapper: $runCmdPath`r`n",
        $utf8NoBom
    )
    $engineLogPath = Join-Path $RunLogDir "icarion-engine-$runStamp.log"
    @(
        "@echo off",
        "setlocal",
        "echo [launcher] wrapper started>> `"$runLogPath`"",
        "cd /d `"$Root`"",
        "echo [launcher] cwd=%CD%>> `"$runLogPath`"",
        "echo [launcher] exe=`"$IcarionExe`">> `"$runLogPath`"",
        "echo [launcher] arg=`"$configArg`">> `"$runLogPath`"",
        "echo [launcher] engine_log=`"$engineLogPath`">> `"$runLogPath`"",
        "echo [launcher] command=`"$IcarionExe`" --log-level DEBUG --log-file `"$engineLogPath`" `"$configArg`">> `"$runLogPath`"",
        "`"$IcarionExe`" --log-level DEBUG --log-file `"$engineLogPath`" `"$configArg`" >> `"$runLogPath`" 2>&1",
        "set `"ICARION_EXIT=%ERRORLEVEL%`"",
        "echo [launcher] exit_code=%ICARION_EXIT%>> `"$runLogPath`"",
        "exit /b %ICARION_EXIT%"
    ) | Set-Content -LiteralPath $runCmdPath -Encoding ASCII
    Append-Log "Generated wrapper:"
    Get-Content -LiteralPath $runCmdPath | ForEach-Object { Append-Log ("  " + $_) }
    Append-Log ""

    $startInfo = New-Object System.Diagnostics.ProcessStartInfo
    $startInfo.FileName = if ($env:ComSpec) { $env:ComSpec } else { "cmd.exe" }
    $startInfo.Arguments = '/d /s /c call "' + $runCmdPath + '"'
    $startInfo.WorkingDirectory = $Root
    $startInfo.UseShellExecute = $false
    $startInfo.CreateNoWindow = $true

    $script:process = New-Object System.Diagnostics.Process
    $script:process.StartInfo = $startInfo

    try {
        Set-RunControlsEnabled $false
        $statusLabel.Text = "Running..."

        [void]$script:process.Start()
        $logTimer.Start()
    } catch {
        $logTimer.Stop()
        Set-RunControlsEnabled $true
        $statusLabel.Text = "Failed to start"
        Append-Log $_.Exception.Message
        Show-Message `
            -Message $_.Exception.Message `
            -Title "Failed to start ICARION" `
            -Icon ([System.Windows.Forms.MessageBoxIcon]::Error)
    }
})

function Start-Analysis {
    param(
        [string]$Name,
        [string]$ScriptName,
        [string[]]$ExtraArgs = @(),
        [string]$CsvOption = "--out-csv",
        [bool]$PreferPerSpeciesPlot = $false,
        [string]$OutputExtension = ".png"
    )

    if (-not (Test-Path $AnalysisDir)) {
        Show-Message `
            -Message "Cannot find packaged analysis scripts.`nExpected:`n$AnalysisDir" `
            -Title "Analysis unavailable" `
            -Icon ([System.Windows.Forms.MessageBoxIcon]::Error)
        return
    }
    $pythonCommand = Get-PythonCommand
    if (-not $pythonCommand) {
        Show-Message `
            -Message "Python was not found. Install Python 3 and the packages from requirements-analysis.txt." `
            -Title "Analysis unavailable" `
            -Icon ([System.Windows.Forms.MessageBoxIcon]::Error)
        return
    }

    $trajPath = Select-TrajectoryFile
    if (-not $trajPath) {
        return
    }

    New-Item -ItemType Directory -Force $RunLogDir | Out-Null
    $runStamp = Get-Date -Format "yyyyMMdd-HHmmss"
    $analysisOutDir = Join-Path $Root "analysis-output"
    New-Item -ItemType Directory -Force $analysisOutDir | Out-Null
    $safeName = $Name.ToLowerInvariant().Replace(" ", "-")
    $runLogPath = Join-Path $RunLogDir "$safeName-$runStamp.log"
    $runCmdPath = Join-Path $RunLogDir "$safeName-$runStamp.cmd"
    $plotPath = Join-Path $analysisOutDir "$safeName-$runStamp$OutputExtension"
    $previewPath = $plotPath
    $perSpeciesPath = Join-Path $analysisOutDir "$safeName-$runStamp-per-species.png"
    if ($PreferPerSpeciesPlot) {
        $previewPath = $perSpeciesPath
    }
    $csvPath = if ($CsvOption) { Join-Path $analysisOutDir "$safeName-$runStamp.csv" } else { "" }
    $scriptPath = Join-Path $AnalysisDir $ScriptName
    $script:runLogPath = $runLogPath
    $script:currentPlotPath = $previewPath
    $script:currentCsvPath = $csvPath
    $script:currentRunName = $Name

    $logBox.Clear()
    Append-Log "Analysis: $Name"
    Append-Log "Trajectory: $trajPath"
    Append-Log "Plot: $plotPath"
    if ($PreferPerSpeciesPlot) {
        Append-Log "Per-species plot: $perSpeciesPath"
    }
    Append-Log "CSV:  $csvPath"
    Append-Log "Log:  $runLogPath"
    Append-Log ""
    $utf8NoBom = New-Object System.Text.UTF8Encoding -ArgumentList $false
    [System.IO.File]::WriteAllText($runLogPath, "Analysis: $Name`r`nTrajectory: $trajPath`r`nPlot: $plotPath`r`nCSV: $csvPath`r`n", $utf8NoBom)

    $extra = ""
    foreach ($arg in $ExtraArgs) {
        $extra += " $arg"
    }
    $csvArgs = ""
    if ($CsvOption) {
        $csvArgs = " $CsvOption `"$csvPath`""
    }
    $plotArgs = ""
    if ($PreferPerSpeciesPlot) {
        $plotArgs = " --out-per-species `"$perSpeciesPath`""
    }
    @(
        "@echo off",
        "setlocal",
        "echo [analysis-launcher] wrapper started>> `"$runLogPath`"",
        "cd /d `"$Root`"",
        "echo [analysis-launcher] cwd=%CD%>> `"$runLogPath`"",
        "set `"PYTHONPATH=$Root;$AnalysisDir;%PYTHONPATH%`"",
        "echo [analysis-launcher] python=$pythonCommand>> `"$runLogPath`"",
        "echo [analysis-launcher] script=`"$scriptPath`">> `"$runLogPath`"",
        "echo [analysis-launcher] trajectory=`"$trajPath`">> `"$runLogPath`"",
        "echo [analysis-launcher] plot=`"$plotPath`">> `"$runLogPath`"",
        "echo [analysis-launcher] preview=`"$previewPath`">> `"$runLogPath`"",
        "echo [analysis-launcher] csv=`"$csvPath`">> `"$runLogPath`"",
        "echo [analysis-launcher] PYTHONPATH=%PYTHONPATH%>> `"$runLogPath`"",
        "echo [analysis-launcher] command=$pythonCommand `"$scriptPath`" --traj `"$trajPath`" --out `"$plotPath`"$csvArgs$plotArgs$extra>> `"$runLogPath`"",
        "$pythonCommand `"$scriptPath`" --traj `"$trajPath`" --out `"$plotPath`"$csvArgs$plotArgs$extra >> `"$runLogPath`" 2>&1",
        "set `"ANALYSIS_EXIT=%ERRORLEVEL%`"",
        "echo [analysis-launcher] exit_code=%ANALYSIS_EXIT%>> `"$runLogPath`"",
        "exit /b %ANALYSIS_EXIT%"
    ) | Set-Content -LiteralPath $runCmdPath -Encoding ASCII
    Append-Log "Generated wrapper:"
    Get-Content -LiteralPath $runCmdPath | ForEach-Object { Append-Log ("  " + $_) }
    Append-Log ""

    $startInfo = New-Object System.Diagnostics.ProcessStartInfo
    $startInfo.FileName = if ($env:ComSpec) { $env:ComSpec } else { "cmd.exe" }
    $startInfo.Arguments = '/d /s /c call "' + $runCmdPath + '"'
    $startInfo.WorkingDirectory = $Root
    $startInfo.UseShellExecute = $false
    $startInfo.CreateNoWindow = $true

    $script:process = New-Object System.Diagnostics.Process
    $script:process.StartInfo = $startInfo

    try {
        Set-RunControlsEnabled $false
        $statusLabel.Text = "Analyzing..."
        [void]$script:process.Start()
        $logTimer.Start()
    } catch {
        $logTimer.Stop()
        Set-RunControlsEnabled $true
        $statusLabel.Text = "Failed to start analysis"
        Append-Log $_.Exception.Message
        Show-Message `
            -Message $_.Exception.Message `
            -Title "Failed to start analysis" `
            -Icon ([System.Windows.Forms.MessageBoxIcon]::Error)
    }
}

$analysisMobilityButton.Add_Click({
    Start-Analysis "IMS mobility" "ims_mobility.py"
})

$analysisArrivalButton.Add_Click({
    Start-Analysis `
        -Name "Arrival times" `
        -ScriptName "arrival_time_distribution.py" `
        -PreferPerSpeciesPlot $true
})

$analysisTrajectoryButton.Add_Click({
    Start-Analysis `
        -Name "Trajectories" `
        -ScriptName "plot_trajectories.py" `
        -ExtraArgs @("--time-stride", "5", "--max-ions", "120", "--max-per-species", "40") `
        -CsvOption ""
})

$analysisMeanButton.Add_Click({
    Start-Analysis `
        -Name "Mean positions" `
        -ScriptName "mean_positions.py" `
        -ExtraArgs @("--time-stride", "5", "--max-ions", "400", "--max-per-species", "200") `
        -CsvOption "--csv"
})

$analysisEliminationButton.Add_Click({
    Start-Analysis `
        -Name "Eliminations" `
        -ScriptName "elimination_histograms.py" `
        -ExtraArgs @("--max-ions", "500", "--max-per-species", "200") `
        -CsvOption ""
})

$analysisAnimateButton.Add_Click({
    Start-Analysis `
        -Name "Animation" `
        -ScriptName "animate_trajectories.py" `
        -ExtraArgs @("--projection", "xy", "--style", "trail", "--theme", "dark", "--max-ions", "80", "--max-per-species", "25", "--time-stride", "4", "--frame-step", "2", "--max-frames", "300", "--writer", "pillow", "--fps", "15", "--dpi", "120") `
        -CsvOption "" `
        -OutputExtension ".gif"
})

$stopButton.Add_Click({
    if ($script:process -and -not $script:process.HasExited) {
        $script:process.Kill()
        $statusLabel.Text = "Stopping..."
    }
})

$form.Add_FormClosing({
    param($sender, $eventArgs)
    if ($script:process -and -not $script:process.HasExited) {
        $result = [System.Windows.Forms.MessageBox]::Show(
            $form,
            "A simulation is still running. Stop it and close?",
            "ICARION Launcher",
            [System.Windows.Forms.MessageBoxButtons]::YesNo,
            [System.Windows.Forms.MessageBoxIcon]::Question
        )
        if ($result -eq [System.Windows.Forms.DialogResult]::Yes) {
            $script:process.Kill()
        } else {
            $eventArgs.Cancel = $true
        }
    }
})

if (Test-Path (Join-Path $DefaultConfigDir "ims\ims_basic.json")) {
    $configBox.Text = Join-Path $DefaultConfigDir "ims\ims_basic.json"
}

[void]$form.ShowDialog()
