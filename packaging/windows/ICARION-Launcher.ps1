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
$configPanel.Height = 88
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
$configBox.Size = New-Object System.Drawing.Size(500, 24)
$configBox.Anchor = "Top,Left,Right"
$configBox.Font = $FontUi
$configBox.BorderStyle = [System.Windows.Forms.BorderStyle]::FixedSingle

$browseButton = New-FlatButton "Browse..." ([System.Drawing.Color]::FromArgb(229, 234, 240)) $ColorText
$browseButton.Location = New-Object System.Drawing.Point(526, 38)
$browseButton.Size = New-Object System.Drawing.Size(86, 30)
$browseButton.Anchor = "Top,Right"

$examplesButton = New-FlatButton "Examples..." ([System.Drawing.Color]::FromArgb(229, 234, 240)) $ColorText
$examplesButton.Location = New-Object System.Drawing.Point(620, 38)
$examplesButton.Size = New-Object System.Drawing.Size(96, 30)
$examplesButton.Anchor = "Top,Right"

$runButton = New-FlatButton "Run" $ColorPrimary
$runButton.Location = New-Object System.Drawing.Point(724, 38)
$runButton.Size = New-Object System.Drawing.Size(68, 30)
$runButton.Anchor = "Top,Right"

$stopButton = New-FlatButton "Stop" $ColorStop
$stopButton.Text = "Stop"
$stopButton.Enabled = $false
$stopButton.Location = New-Object System.Drawing.Point(800, 38)
$stopButton.Size = New-Object System.Drawing.Size(60, 30)
$stopButton.Anchor = "Top,Right"

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

$contentPanel.Controls.Add($logPanel)
$contentPanel.Controls.Add($configPanel)

$form.Controls.Add($contentPanel)
$form.Controls.Add($statusPanel)
$form.Controls.Add($headerPanel)

$script:process = $null

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

    $logBox.Clear()
    Append-Log "ICARION: $IcarionExe"
    Append-Log "Config:  $configPath"
    Append-Log "Log:     $runLogPath"
    Append-Log ""
    "ICARION: $IcarionExe`r`nConfig:  $configPath`r`n" | Set-Content -LiteralPath $runLogPath

    $startInfo = New-Object System.Diagnostics.ProcessStartInfo
    $startInfo.FileName = $IcarionExe
    $startInfo.Arguments = '"' + $configPath + '"'
    $startInfo.WorkingDirectory = $Root
    $startInfo.UseShellExecute = $false
    $startInfo.RedirectStandardOutput = $true
    $startInfo.RedirectStandardError = $true
    $startInfo.CreateNoWindow = $true

    $script:process = New-Object System.Diagnostics.Process
    $script:process.StartInfo = $startInfo
    $script:process.EnableRaisingEvents = $true

    $outputHandler = [System.Diagnostics.DataReceivedEventHandler]{
        param($sender, $eventArgs)
        if ($eventArgs.Data -ne $null) {
            $line = $eventArgs.Data
            Add-Content -LiteralPath $runLogPath -Value $line
            $form.BeginInvoke([Action]{ Append-Log $line }) | Out-Null
        }
    }
    $errorHandler = [System.Diagnostics.DataReceivedEventHandler]{
        param($sender, $eventArgs)
        if ($eventArgs.Data -ne $null) {
            $line = $eventArgs.Data
            Add-Content -LiteralPath $runLogPath -Value $line
            $form.BeginInvoke([Action]{ Append-Log $line }) | Out-Null
        }
    }
    $exitHandler = [System.EventHandler]{
        param($sender, $eventArgs)
        $code = $sender.ExitCode
        $form.BeginInvoke([Action]{
            $runButton.Enabled = $true
            $browseButton.Enabled = $true
            $examplesButton.Enabled = $true
            $stopButton.Enabled = $false
            $statusLabel.Text = "Finished with exit code $code"
            Append-Log ""
            Append-Log "Finished with exit code $code"
            if ($code -ne 0) {
                Show-Message `
                    -Message "ICARION stopped with exit code $code.`n`nThe log was saved here:`n$runLogPath" `
                    -Title "Run failed" `
                    -Icon ([System.Windows.Forms.MessageBoxIcon]::Error)
            }
        }) | Out-Null
    }

    $script:process.add_OutputDataReceived($outputHandler)
    $script:process.add_ErrorDataReceived($errorHandler)
    $script:process.add_Exited($exitHandler)

    try {
        $runButton.Enabled = $false
        $browseButton.Enabled = $false
        $examplesButton.Enabled = $false
        $stopButton.Enabled = $true
        $statusLabel.Text = "Running..."

        [void]$script:process.Start()
        $script:process.BeginOutputReadLine()
        $script:process.BeginErrorReadLine()
    } catch {
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
