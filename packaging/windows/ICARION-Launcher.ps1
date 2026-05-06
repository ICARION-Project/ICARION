Add-Type -AssemblyName System.Windows.Forms
Add-Type -AssemblyName System.Drawing

$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$IcarionExe = Join-Path $Root "bin\icarion.exe"
$DefaultConfigDir = Join-Path $Root "examples"

function Append-Log {
    param([string]$Text)
    if ([string]::IsNullOrWhiteSpace($Text)) {
        return
    }
    $logBox.AppendText($Text + [Environment]::NewLine)
    $logBox.SelectionStart = $logBox.Text.Length
    $logBox.ScrollToCaret()
}

$form = New-Object System.Windows.Forms.Form
$form.Text = "ICARION Launcher"
$form.StartPosition = "CenterScreen"
$form.Size = New-Object System.Drawing.Size(900, 620)
$form.MinimumSize = New-Object System.Drawing.Size(720, 460)

$configLabel = New-Object System.Windows.Forms.Label
$configLabel.Text = "Config file"
$configLabel.Location = New-Object System.Drawing.Point(12, 18)
$configLabel.AutoSize = $true

$configBox = New-Object System.Windows.Forms.TextBox
$configBox.Location = New-Object System.Drawing.Point(90, 14)
$configBox.Size = New-Object System.Drawing.Size(610, 24)
$configBox.Anchor = "Top,Left,Right"

$browseButton = New-Object System.Windows.Forms.Button
$browseButton.Text = "Browse..."
$browseButton.Location = New-Object System.Drawing.Point(710, 12)
$browseButton.Size = New-Object System.Drawing.Size(80, 28)
$browseButton.Anchor = "Top,Right"

$runButton = New-Object System.Windows.Forms.Button
$runButton.Text = "Run"
$runButton.Location = New-Object System.Drawing.Point(798, 12)
$runButton.Size = New-Object System.Drawing.Size(74, 28)
$runButton.Anchor = "Top,Right"

$stopButton = New-Object System.Windows.Forms.Button
$stopButton.Text = "Stop"
$stopButton.Enabled = $false
$stopButton.Location = New-Object System.Drawing.Point(798, 46)
$stopButton.Size = New-Object System.Drawing.Size(74, 28)
$stopButton.Anchor = "Top,Right"

$statusLabel = New-Object System.Windows.Forms.Label
$statusLabel.Text = "Ready"
$statusLabel.Location = New-Object System.Drawing.Point(12, 52)
$statusLabel.Size = New-Object System.Drawing.Size(760, 22)
$statusLabel.Anchor = "Top,Left,Right"

$logBox = New-Object System.Windows.Forms.TextBox
$logBox.Location = New-Object System.Drawing.Point(12, 86)
$logBox.Size = New-Object System.Drawing.Size(860, 480)
$logBox.Anchor = "Top,Bottom,Left,Right"
$logBox.Multiline = $true
$logBox.ScrollBars = "Both"
$logBox.WordWrap = $false
$logBox.ReadOnly = $true
$logBox.Font = New-Object System.Drawing.Font("Consolas", 9)

$form.Controls.AddRange(@(
    $configLabel,
    $configBox,
    $browseButton,
    $runButton,
    $stopButton,
    $statusLabel,
    $logBox
))

$script:process = $null

$browseButton.Add_Click({
    $dialog = New-Object System.Windows.Forms.OpenFileDialog
    $dialog.Title = "Select ICARION config"
    $dialog.Filter = "JSON config (*.json)|*.json|All files (*.*)|*.*"
    if (Test-Path $DefaultConfigDir) {
        $dialog.InitialDirectory = $DefaultConfigDir
    }
    if ($dialog.ShowDialog($form) -eq [System.Windows.Forms.DialogResult]::OK) {
        $configBox.Text = $dialog.FileName
    }
})

$runButton.Add_Click({
    if (-not (Test-Path $IcarionExe)) {
        [System.Windows.Forms.MessageBox]::Show(
            $form,
            "Cannot find bin\icarion.exe next to this launcher.",
            "ICARION Launcher",
            [System.Windows.Forms.MessageBoxButtons]::OK,
            [System.Windows.Forms.MessageBoxIcon]::Error
        ) | Out-Null
        return
    }

    $configPath = $configBox.Text.Trim()
    if (-not (Test-Path $configPath)) {
        [System.Windows.Forms.MessageBox]::Show(
            $form,
            "Please select a valid JSON config file.",
            "ICARION Launcher",
            [System.Windows.Forms.MessageBoxButtons]::OK,
            [System.Windows.Forms.MessageBoxIcon]::Warning
        ) | Out-Null
        return
    }

    $logBox.Clear()
    Append-Log "ICARION: $IcarionExe"
    Append-Log "Config:  $configPath"
    Append-Log ""

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
            $form.BeginInvoke([Action]{ Append-Log $line }) | Out-Null
        }
    }
    $errorHandler = [System.Diagnostics.DataReceivedEventHandler]{
        param($sender, $eventArgs)
        if ($eventArgs.Data -ne $null) {
            $line = $eventArgs.Data
            $form.BeginInvoke([Action]{ Append-Log $line }) | Out-Null
        }
    }
    $exitHandler = [System.EventHandler]{
        param($sender, $eventArgs)
        $code = $sender.ExitCode
        $form.BeginInvoke([Action]{
            $runButton.Enabled = $true
            $browseButton.Enabled = $true
            $stopButton.Enabled = $false
            $statusLabel.Text = "Finished with exit code $code"
            Append-Log ""
            Append-Log "Finished with exit code $code"
        }) | Out-Null
    }

    $script:process.add_OutputDataReceived($outputHandler)
    $script:process.add_ErrorDataReceived($errorHandler)
    $script:process.add_Exited($exitHandler)

    try {
        $runButton.Enabled = $false
        $browseButton.Enabled = $false
        $stopButton.Enabled = $true
        $statusLabel.Text = "Running..."

        [void]$script:process.Start()
        $script:process.BeginOutputReadLine()
        $script:process.BeginErrorReadLine()
    } catch {
        $runButton.Enabled = $true
        $browseButton.Enabled = $true
        $stopButton.Enabled = $false
        $statusLabel.Text = "Failed to start"
        Append-Log $_.Exception.Message
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
