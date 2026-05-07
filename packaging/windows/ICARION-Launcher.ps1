Add-Type -AssemblyName System.Windows.Forms
Add-Type -AssemblyName System.Drawing

$ErrorActionPreference = "Stop"
[System.Windows.Forms.Application]::EnableVisualStyles()
[System.Windows.Forms.Application]::SetCompatibleTextRenderingDefault($false)

$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$IcarionExe = Join-Path $Root "bin\icarion.exe"
$DefaultConfigDir = Join-Path $Root "examples"
$RunLogDir = Join-Path $Root "launcher-logs"

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
$configPanel.Height = 126
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

$configPanel.Controls.AddRange(@(
    $configLabel,
    $configBox,
    $browseButton,
    $examplesButton,
    $runButton,
    $stopButton
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
            $runButton.Enabled = $true
            $browseButton.Enabled = $true
            $examplesButton.Enabled = $true
            $stopButton.Enabled = $false
            $statusLabel.Text = "Finished with exit code $code"
            Append-Log ""
            Append-Log "Finished with exit code $code"
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
        $runButton.Enabled = $true
        $browseButton.Enabled = $true
        $examplesButton.Enabled = $true
        $stopButton.Enabled = $false
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
    $script:runLogPath = $runLogPath

    $logBox.Clear()
    Append-Log "ICARION: $IcarionExe"
    Append-Log "Config:  $configPath"
    Append-Log "Log:     $runLogPath"
    Append-Log ""
    $utf8NoBom = New-Object System.Text.UTF8Encoding -ArgumentList $false
    [System.IO.File]::WriteAllText($runLogPath, "ICARION: $IcarionExe`r`nConfig:  $configPath`r`n", $utf8NoBom)
    @(
        "@echo off",
        "cd /d `"$Root`"",
        "`"$IcarionExe`" `"$configPath`" >> `"$runLogPath`" 2>&1",
        "exit /b %ERRORLEVEL%"
    ) | Set-Content -LiteralPath $runCmdPath -Encoding ASCII

    $startInfo = New-Object System.Diagnostics.ProcessStartInfo
    $startInfo.FileName = if ($env:ComSpec) { $env:ComSpec } else { "cmd.exe" }
    $startInfo.Arguments = '/d /c "' + $runCmdPath + '"'
    $startInfo.WorkingDirectory = $Root
    $startInfo.UseShellExecute = $false
    $startInfo.CreateNoWindow = $true

    $script:process = New-Object System.Diagnostics.Process
    $script:process.StartInfo = $startInfo

    try {
        $runButton.Enabled = $false
        $browseButton.Enabled = $false
        $examplesButton.Enabled = $false
        $stopButton.Enabled = $true
        $statusLabel.Text = "Running..."

        [void]$script:process.Start()
        $logTimer.Start()
    } catch {
        $logTimer.Stop()
        $runButton.Enabled = $true
        $browseButton.Enabled = $true
        $examplesButton.Enabled = $true
        $stopButton.Enabled = $false
        $statusLabel.Text = "Failed to start"
        Append-Log $_.Exception.Message
        Show-Message `
            -Message $_.Exception.Message `
            -Title "Failed to start ICARION" `
            -Icon ([System.Windows.Forms.MessageBoxIcon]::Error)
    }
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
